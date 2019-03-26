function [activity_clusters] = clustering_k_means(calcium_data,struct,varargin)
%% February 2019 - Carey lab - leonard.dupont@ens.fr
%..........................................................................
% This whole script is dedicated to an attempt of spatially clustering the
% cells based on their calcium activity.
% ........................................................................
% 
% -- INPUT ------------
%
%    struct            struct with roi data
%
%    calcium_data      a n x p matrix where (n) are the rois and (p) the
%                      observations (time: fluorescence value for the roi)
%
%
%    varargin (opt.)   K            nb of clusters for the kmean clustering
%                      runs         nb of times that the kmeans will run
%                      P            threshold probability of occurence to 
%                                   keep roi in neighbourhood
%                      Dth          threshold distance above with clusters 
%                                   are not merged anymore: 
%                                   (%common elements)^{-1}
%                      spatialplot  1 if you want it, 0 otherwise (default)
%                      rasterplot   1 if you want it, 0 otherwise (default)
%                      frame_path   the path to the registration template
%                                   (if spatialplot == 1)
%                      deconvolutO  1 if you use the deconvolution struct
%
% -- OUTPUT -----------
%
%    activity_clusters      struct with fields
%                           >    cluster_i    >  roi_j > centroid [x,y]
%                                  (...)               > spiketimes
%                                                      > raw_trace (Fluo.)                        
%                                             >  (...)
%                         
% ........................................................................



ip = inputParser;

ip.addParameter('K',5);
ip.addParameter('runs',2500);
ip.addParameter('P',0.8);
ip.addParameter('Dth',1/0.6);
ip.addParameter('spatialplot',0);
ip.addParameter('rasterplot',0);
ip.addParameter('frame_path',[]);
ip.addParameter('deconv0',0);

parse(ip,varargin{:});

K = ip.Results.K;
runs = ip.Results.runs;
P = ip.Results.P;
Dth = ip.Results.Dth;
spatialplot = logical(ip.Results.spatialplot);
rasterplot = logical(ip.Results.rasterplot);
frame_path = ip.Results.frame_path;
deconv0 = ip.Results.deconv0;

%% normalise and bring back to baseline (0)
tic
if deconv0
    u = numel(fieldnames(struct));
else
    u = struct.n_cells;
end
disp(' ------ 1 - Normalising the calcium data ------')
%zm_calcium_data = zero_and_mean(calcium_data); %normalise with the mean
zm_calcium_data = zero_and_max(calcium_data);  %normalise with the max 
zm_calcium_data = zm_calcium_data(1:u,:);
[rois,datapoints] = size(zm_calcium_data);
toc
%% K-MEANS ++ algorithm
clear clusterresults
disp([' ------ 2 - Running the k-means++ algorithm ',num2str(runs),' times with ',num2str(K),' centroids ------'])

for i = 1:runs 
    idx = kmeans(zm_calcium_data, K,'Distance','correlation'); %running k-means with
    for N = 1:K
        clusterresults.(['run_',num2str(i)]).(['c',num2str(N)]) = find(idx==N); 
    end
end 
toc
%% Neighbouring
% We now have K clusters replicated 2500 times. We want to see how often an
% ROI has been clustered with others, so we will first build a big vector
% for each ROI containing the cells it's been put with. 
%%
disp(' ------ 3 - Finding neighbours ------')
clear frequencies

% Hereunder, we take the first run and build a struct called 'frequencies'
% with subfields being individual rois. In each ROIs' 'neighbours' vector,
% we store the other ROIs it was with (i.e the cluster c containing the ROI
% in run 1).  

run = 1;
first_clusters = clusterresults.(['run_',num2str(run)]);

for N = 1:K
    cluster = first_clusters.(['c',num2str(N)]);
    for roi = 1:length(cluster)
        region = cluster(roi);
        frequencies.(['roi_',num2str(region)]).neighbours = cluster;
    end
end

%%
% Now for each of the remaining runs (from 2 to Nruns), we find in which
% cluster each ROI was and store the neighbours again (i.e the cluster c
% containing the ROI in run i). 

for i = 2:runs
    for roi = 1:rois
        c = 1;
        cluster = clusterresults.(['run_',num2str(i)]).(['c',num2str(c)]);
        while isempty(find(cluster == roi))
            c = c + 1;
            cluster = clusterresults.(['run_',num2str(i)]).(['c',num2str(c)]);
        end
            konkat = cat(1,frequencies.(['roi_',num2str(roi)]).neighbours, cluster);
            frequencies.(['roi_',num2str(roi)]).neighbours = konkat;    
    end
end
toc
%%
% At this point, we have a huge concatenation vector for each of our ROIs
% in which we can find the ROIs it was sorted with during the k-means. What
% we now want is to prune this list to a stable ensemble : the
% neighbourhood. Only ROIs that have been neighbours in more than 80% of
% the run will be accepted in the neighbourhood (it's an elitist circle). 
Z = 10;

disp([' ------ 4 - Reducing to neighbourhoods with ',num2str(P),' occurence probability ------'])
clear neighbourhood


for roi = 1:rois
    neigh = frequencies.(['roi_',num2str(roi)]).neighbours;
    neighB = [];
    l = length(neigh);
    for k = 1:rois
        probability = length( find(neigh == k) ) / runs;
        if probability > P
            neighB(end+1) = k;
        end
        if length(neighB) > Z
            neighbourhood.(['roi_',num2str(roi)]).neighB = neighB;
        else
            neighbourhood.(['roi_',num2str(roi)]).neighB = roi;
        end
    end
end
toc
%%
% We have now started building a new struct called 'neighbourhood'.
% Subfields are ROI names once again. For each ROI, we have a vector
% 'neighB' comprising the above-selected neighbourhood. 
% As explained in the pdf, we also need a weighing vector for each
% neighbourhood to calculate distances. We create the original weights
% here:
disp(' ------ 5 - Building weight vectors ------')
for roi = 1:rois
    weights = zeros(rois,1);
    itsneighbours = neighbourhood.(['roi_',num2str(roi)]).neighB;
    for k = 1:rois
        if ~isempty( find(itsneighbours == k) ) 
            weights(k) = 1;
        end
    end
    neighbourhood.(['roi_',num2str(roi)]).weights = weights;
end
toc
%%
% We now have to fields per ROI in our neighbourhood struct:
%
%  neighbourhood    >   roi_1    >    neighB : stable neighbours
%                                >    W = [w1, ... , wR], wi = 1 if 
%                                roi i belongs to neighB, 0 otherwise (for
%                                now)
%                   >   roi_2    >>
%                     (...) 
%                   >   roi_R    >>
%
% Hereunder, we start calculating distance matrices and fusing the close
% neighB vectors to form bigger neighbourhoods. The maths are explained in
% the pdf file. 

Dmin = 1;
merge_events = 0;

disp([' ------ 6 - Now merging using hierarchical clustering, parameters Dmin = ',num2str(Dmin),' and Dthresh = ',num2str(Dth),'------'])

while Dmin ~= 987654321
    
    L = numel(fieldnames(neighbourhood));
    names = fieldnames(neighbourhood);
    M = zeros(L,L);
   
    for i = 1:L
        for j = 1:L
            
            if i~=j
                R1 = neighbourhood.(names{i}).neighB;
                R2 = neighbourhood.(names{j}).neighB;

                w1 = neighbourhood.(names{i}).weights;
                w2 = neighbourhood.(names{j}).weights;

                M(i,j) = distance_star(R1,R2,w1,w2);
                
            else
                M(i,j) = 987654321;
            end

        end
    end

    [Dmin,neighbourhood,merge_events] = merge_closest(M,neighbourhood,rois, merge_events,Dth); 
end
toc
%%
disp(' ------ 7 - Finding redundant elements and sorting them based on weights -------')
L = numel(fieldnames(neighbourhood));
names = fieldnames(neighbourhood);

for k = 1:L
    for j = 1:L
         if k~=j
             R1 = neighbourhood.(names{k}).neighB;
             R2 = neighbourhood.(names{j}).neighB;

             w1 = neighbourhood.(names{k}).weights;
             w2 = neighbourhood.(names{j}).weights;

             R1nR2 = intersect(R1,R2);
             inter = length(R1nR2);

             if inter > 0
                 for i = 1:inter
                     redundant = R1nR2(i);
                     if w1(redundant) ~= w2(redundant)
                         closer = max(w1(redundant),w2(redundant));
                         if closer == w1(redundant)
                             index = find(R2 == redundant);
                             R2(index) = [];
                             w2(index) = 0;
                         elseif closer == w2(redundant)
                             index = find(R1 == redundant);
                             R1(index) = [];
                             w1(index) = 0;
                         end
                     end
                 end
             end

             neighbourhood.(names{k}).neighB = R1;
             neighbourhood.(names{j}).neighB = R2;

             neighbourhood.(names{k}).weights = w1;
             neighbourhood.(names{k}).weights = w2;
         end
    end
end
toc
%% Removing the clusters with fewer than Z cells
L = numel(fieldnames(neighbourhood));
names = fieldnames(neighbourhood);
Z = 5;
disp([' ------ 8 - Now removing clusters comprising less than ',num2str(Z),' cells --------'])
disp(['_i_ Initial situation = ',num2str(L),' clusters.'])

for k = 1:L
    x = length(neighbourhood.(names{k}).neighB);
    if x < Z
        neighbourhood = rmfield(neighbourhood,names{k});
    end
end
L = numel(fieldnames(neighbourhood));
disp(['_ii_ Final situation = ',num2str(L),' clusters'])

%% Now we build a final struct with the clusters but with activity data AND spatial coordinates (centroids)
disp(' ------ 9 - Preparing spatial plot, rasterplot and output structure --------')
clear activity_clusters
L = numel(fieldnames(neighbourhood));
names = fieldnames(neighbourhood);

for k = 1:L
    regions = neighbourhood.(names{k}).neighB;
    H = length(regions);
    x = zeros(H,1);
    y = zeros(H,1);
    calcium_matrix = zeros(H,datapoints);
    
    if deconv0
        deconvolution = struct;
        for roi = 1:H
            nb = regions(roi);
            activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(nb)]).centroid = deconvolution.(['roi_',num2str(nb)]).centroid;
            x(roi) = deconvolution.(['roi_',num2str(nb)]).centroid(1,1);
            y(roi) = deconvolution.(['roi_',num2str(nb)]).centroid(1,2);
            activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(nb)]).spiketimes = deconvolution.(['roi_',num2str(nb)]).spiketimes;
            activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(nb)]).raw_trace = deconvolution.(['roi_',num2str(nb)]).raw_trace;
            calcium_matrix(roi,:) = (deconvolution.(['roi_',num2str(nb)]).raw_trace).';
        end
        
    else
        cn = struct;
        for roi = 1:H
            nb = regions(roi);
            activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(nb)]).centroid = cn.centroid{1,nb};
            x(roi) = activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(nb)]).centroid(1,1);
            y(roi) = activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(nb)]).centroid(1,2);
            activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(nb)]).raw_trace = cn.intensity(:,nb);
            calcium_matrix(roi,:) = cn.intensity(:,nb);
        end
    end
    
    activity_clusters.(['cluster_',num2str(k)]).mean_trace = mean(calcium_matrix,1);
    xc = sum(x) / H;
    yc = sum(y) / H;
    activity_clusters.(['cluster_',num2str(k)]).barycenter(1,1:2) = [xc yc];
end

%% Correlation matrix

C = zeros(L,L);

for k = 1:L
    for j = 1:L
        Mt1 = activity_clusters.(['cluster_',num2str(k)]).mean_trace;
        Mt2 = activity_clusters.(['cluster_',num2str(j)]).mean_trace;
        X = corrcoef(Mt1,Mt2); 
        C(k,j) = X(1,2);
    end
end

%% We can now plot the rois back on a frame with a colour-based code in compliance with their cluster
if spatialplot  
    scatter_clusters(activity_clusters,frame_path) 
end

figure; hold on
imagesc(C), colorbar
hold off 

%% We can also give back a raster plot
if rasterplot
   raster_clusters(activity_clusters); 
   %only works if you fed the clustering algorithm with deconvolution
end
end






