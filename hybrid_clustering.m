function [activity_clusters,C] = hybrid_clustering(data,cl)
%% March 2019 - Carey lab - leonard.dupont@ens.fr
%..........................................................................
% This whole script is dedicated to an attempt of spatially clustering the
% cells based on their calcium activity. It can work on either cn neuron
% structure or on tif movies, in which case it will use a pixel-based 
% clustering method. 
% ........................................................................
% 
% -- INPUT ------------
%
%    data   (i)   either cn struct with AT LEAST following fields
%                       n_cells (number of cells)
%                       intensity (usually matrix(timeframes,cells))
%                       centroid (cell)
%                       spiketimes if you want to use the rasterplot
%           (ii)  or a tiff movie of the FOV. You can give the path or
%                 the (h * w * time) matrix directly if it is in your wd.
%
%    cl     cl(usterparameters ) struct with fields
%                   
%               K        number of clusters for the kmeans  (default : 5)
%               runs     number of kmeanr runs (default : 2500)
%               p        redundance threshold to create neighbourhoods
%                        (default : 0.8)
%               Dth      threshold distance to merge neighbourhoods
%                        (default : 1/0.6)
%        spatialplot     logical, 1 if you want to plot (default : 0)
%         rasterplot     logical, 1 if you want to plot (default : 0)
%          normalise     logical, 1 if you want to normalise (default : 0)
%        norm_method     either 'mean' or 'max' (default : [])
%               Z        minimum cluster size to be kept (discarded otherw)
%                        (default : 10)
%           gridbins     when using pixel-based clustering, number of bins
%                        to create in x and y directions (default : 15)
%
%
% -- OUTPUT -----------
%
%    activity_clusters      struct with fields
%                           >    cluster_i    >  roi_j > centroid [x,y]
%                                  (...)               > spiketimes
%                                                      > raw_trace (Fluo.)                        
%                                             >  (...)
%                   
%                    C      correlation matrix               

%% function parameters 

cl = check_empty(cl,'K',5);
cl = check_empty(cl,'runs',2500);
cl = check_empty(cl,'p',0.8);
cl = check_empty(cl,'Dth',1/0.6);
cl = check_empty(cl,'spatialplot',0);
cl = check_empty(cl,'rasterplot',0);
cl = check_empty(cl,'normalise',0);
cl = check_empty(cl,'norm_method','max');
cl = check_empty(cl,'Z',10);
cl = check_empty(cl,'gridbins',15);


if ~isstruct(data)
    disp('Using pixels as clustering units.')
    grdfov = grid_grouping(data,cl.gridbins);
    N = (cl.gridbins)^2;
    pixels = 1;
      
elseif isstruct(data)
    disp('Using rois as clustering units.')
    cn = data;
    N = cn.n_cells;
    pixels = 0;
else
    error('Data and method do not match : check documentation. Aborting.')
end



%% Clustering units : regions of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~pixels 
    %% 1 - Normalisation 

    tic

    if logical(cl.normalise)
        disp(' ------ 1 - Normalising the calcium data ------')
        switch cl.norm_method
            case 'mean'
                zm_calcium_data = zero_and_mean(cn.intensity.'); 
            case 'max'
                zm_calcium_data = zero_and_max(cn.intensity.');
        end 
    else
        zm_calcium_data = cn.intensity.'; 
    end
  
%% Clustering units : pixels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif pixels
    %% 1 - Normalisation 

    tic

    if logical(cl.normalise)
        disp(' ------ 1 - Normalising the calcium data ------')
        switch cl.norm_method
            case 'mean'
                zm_calcium_data = zero_and_mean(grdfov.intensity); 
            case 'max'
                zm_calcium_data = zero_and_max(grdfov.intensity);
        end 
    else
        zm_calcium_data = grdfov.intensity; 
    end
    
end  
    %% 2 - Repeated, pseudo-random K-means

    clear clusterresults
    disp([' ------ 2 - Running the k-means++ algorithm ',num2str(cl.runs),...
        ' times with ',num2str(cl.K),' centroids ------'])

    for r = 1:cl.runs 
        idx = kmeans(zm_calcium_data, cl.K,'Distance','correlation'); 
        for k = 1:cl.K
            clusterresults.(['run_',num2str(r)]).(['c',num2str(k)]) = ...
                                                              find(idx==k); 
        end
    end 

    toc

    % clusterresults    >   run_1    >   c1    > roi numbers of cells in c1
    %                       (...)       (...)
    %                                    cK    > roi nbs of cells in cK           
    %                   > run_(runs) 
      

    %% 3 - Neighbouring
    % We now have K clusters replicated 2500 times. We want to see how 
    %often an ROI has been clustered with others, so we will first build 
    % a big vector for each ROI containing the cells it's been put with. 
    
    disp(' ------ 3 - Finding neighbours ------')
    clear frequencies

    % Hereunder, we take the first run and build a struct called 
    %'frequencies' with subfields being individual rois. In each ROIs' 
    %'neighbours' vector, we store the other ROIs it was with 
    % (i.e the cluster c containing the ROI in run 1).  

    run = 1;
    first_clusters = clusterresults.(['run_',num2str(run)]);

    for k = 1:cl.K
        cluster = first_clusters.(['c',num2str(k)]);
        for roi = 1:length(cluster)
            region = cluster(roi);
            frequencies.(['roi_',num2str(region)]).neighbours = cluster;
        end
    end    
    
    % Now for each of the remaining runs (from 2 to Nruns), we find in 
    % which cluster each ROI was and store the neighbours again 
    % (i.e the cluster c containing the ROI in run i). 

    for i = 2:cl.runs 
        for roi = 1:N
            c = 1;
            cluster = ...
                clusterresults.(['run_',num2str(i)]).(['c',num2str(c)]);
            while isempty(find(cluster == roi))
                c = c + 1;
                cluster = ...
                  clusterresults.(['run_',num2str(i)]).(['c',num2str(c)]);
            end
            konkat = ...
            cat(1,frequencies.(['roi_',num2str(roi)]).neighbours, cluster);
            frequencies.(['roi_',num2str(roi)]).neighbours = konkat;    
        end
    end
    toc

    
    %   frequencies  >  roi_1   >  neighbours  = [every roi it was in
    %                   (...)                     during the runs]
    %                >  roi_N                                               
    
   %% 4 - Building 80% confidence neighbourhoods
   % At this point, we have a huge concatenation vector for each of our 
   % ROIs in which we can find the ROIs it was sorted with during the 
   % k-means. What we now want is to prune this list to a stable ensemble :
   % the neighbourhood. Only ROIs that have been neighbours in more than 
   % 80% of the run will be accepted in the neighbourhood 
   % (it is an elitist circle). 
    clear neighbourhood

    disp([' ------ 4 - Reducing to neighbourhoods with ',num2str(cl.p),...
                                          ' occurence probability ------'])
    


    for roi = 1:N
        neighbours = frequencies.(['roi_',num2str(roi)]).neighbours;
        neighB = [];
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        for k = 1:N
            probability = length( find(neighbours == k) ) / cl.runs;
            if probability > cl.p
                neighB(end+1) = k;
            end 
        end
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   
        if length(neighB) > cl.Z
           neighbourhood.(['roi_',num2str(roi)]).neighB = neighB;
        else
           neighbourhood.(['roi_',num2str(roi)]).neighB = roi;
        end 
    end 
    
    %% 5 - Building weight vectors
    % We have now started building a new struct called 'neighbourhood'.
    % Subfields are ROI names once again. For each ROI, we have a vector
    % 'neighB' comprising the above-selected neighbourhood. 
    % As explained in the pdf, we also need a weighing vector for each
    % neighbourhood to calculate distances. We create the original weights
    % here:
    disp(' ------ 5 - Building weight vectors ------')
    
    for roi = 1:N
        weights = zeros(N,1);
        itsneighbours = neighbourhood.(['roi_',num2str(roi)]).neighB;
        for k = 1:rois
            if ~isempty( find(itsneighbours == k) ) 
                weights(k) = 1;
            end
        end
        neighbourhood.(['roi_',num2str(roi)]).weights = weights;
    end
    toc
    
    %% 6 - Hierarchy-based and weighed neighbourhood grouping
    % We now have to fields per ROI in our neighbourhood struct:
    %
    %  neighbourhood    >   roi_1    >    neighB : stable neighbours
    %                                >    W = [w1, ... , wR], wi = 1 if 
    %                                 roi i belongs to neighB, 0 otherwise 
    %                     (...) 
    %                   >   roi_R    
    %
    % Hereunder, we start calculating distance matrices and fusing the 
    % close neighB vectors to form bigger neighbourhoods. 


    Dmin = 1;
    merge_events = 0;

    disp([' ------ 6 - Now merging using hierarchical clustering, parameters Dmin = ',num2str(Dmin),' and Dthresh = ',num2str(cl.Dth),'------'])

    while Dmin ~= 987654321
    % while there exists a distance smaller than Dth...
        
        L = numel(fieldnames(neighbourhood));
        names = fieldnames(neighbourhood);
        M = zeros(L,L);

        for i = 1:L
            for j = 1:L

                if i~=j %not on the diagonal
                    R1 = neighbourhood.(names{i}).neighB;
                    R2 = neighbourhood.(names{j}).neighB;

                    w1 = neighbourhood.(names{i}).weights;
                    w2 = neighbourhood.(names{j}).weights;

                    M(i,j) = distance_star(R1,R2,w1,w2);

                else
                    M(i,j) = 987654321; % M(i,i) = \infty forall i
                end

            end
        end
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        [Dmin,neighbourhood,merge_events] = ...
        merge_closest(M,neighbourhood,rois, merge_events,cl.Dth); 
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    end
    toc    
    
    %% 7 - Some ROIs are redundant between neighBs
    % Hence we remove them from the neighbourhood in which they have the 
    % lowest weight (in which they were the least redundant during merging)
  
    
    disp(' ------ 7 - Finding redundant elements and sorting them based on weights -------')
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    neighbourhood = neighB_redundancy_solver(neighbourhood);
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    toc
    
    
    
    %% 8 - Removing abnormally small clusters
    % Now, a lot of clusters only comprise one cell or two. These 
    % would represent artefacts for most and we ought to deal with it. We
    % can either assign them to the closest neighbourhood (mean
    % fluorescence trace as a distance proxy) OR put them in a special
    % list. 
   

    disp([' ------ 8 - Now removing clusters comprising less than ',num2str(cl.Z),' cells --------'])
    
    L = numel(fieldnames(neighbourhood));
    names = fieldnames(neighbourhood);
    disp(['i. Initial situation = ',num2str(L),' clusters.'])
    
    neighbourhood.lonely = [];
    
    for k = 1:L
        d = length(neighbourhood.(names{k}).neighB);
        if d < cl.Z
            %konkat = ... THINK ABOUT THIS
            %cat(1,neighbourhood.lonely,neighbourhood.(names{k}).neighB);
            %neighbourhood.lonely = konkat;
            neighbourhood = rmfield(neighbourhood,names{k});
        end
    end
    
    L = numel(fieldnames(neighbourhood)); 
    disp(['ii. Final situation = ',num2str(L),' clusters.'])
   
    %% 9 - Building output structure
    disp(' ------ 9 - Preparing output structure and correlation matrix -------')
%% Clustering units : rois %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ~pixels
    clear activity_clusters
    L = numel(fieldnames(neighbourhood));
    names = fieldnames(neighbourhood);

    for k = 1:L
        
        neighB = neighbourhood.(names{k}).neighB;
        neighB_l = length(neighB);
        
        %for barycenter 
        x = zeros(neighB_l,1);
        y = zeros(neighB_l,1);
        
        %for mean calcium trace
        cluster_traces = zeros(neighB_l,datapoints);
        
        for e = 1:neighB_l
            
            roi = neighB(e);
            
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            position = cn.centroid{1,roi};
            activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(roi)]).centroid = position;
            x(e) = position(1,1);
            y(e) = position(1,2);
            clear position
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            ca_data = cn.intensity(:,roi);
            activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(roi)]).intensity = ca_data;
            cluster_traces(e,:) = ca_data;
            clear ca_data
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
            if logical(cl.deconv0)
                spk = cn.spiketimes{1,roi};
                activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(roi)]).spiketimes = spk;
                clear spk
            end
        end
        
        
        activity_clusters.(['cluster_',num2str(k)]).mean_trace = mean(cluster_traces,1);
        xc = sum(x) / neighB_l;
        yc = sum(y) / neighB_l;
        activity_clusters.(['cluster_',num2str(k)]).barycenter(1,1:2) = [xc yc];
  
    end
%% Clustering units : pixels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif pixels
    clear activity_clusters
    L = numel(fieldnames(neighbourhood));
    names = fieldnames(neighbourhood);

    for k = 1:L
        
        neighB = neighbourhood.(names{k}).neighB;
        neighB_l = length(neighB);
        
        %for barycenter 
        x = zeros(neighB_l,1);
        y = zeros(neighB_l,1);
        
        %for mean calcium trace
        cluster_traces = zeros(neighB_l,datapoints);
        
        for e = 1:neighB_l
            
            roi = neighB(e);
            
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            position = grdfov.centroid{1,roi};
            activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(roi)]).centroid = position;
            x(e) = position(1,1);
            y(e) = position(1,2);
            clear position
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            ca_data = grdfov.intensity(:,roi);
            activity_clusters.(['cluster_',num2str(k)]).rois.(['roi_',num2str(roi)]).intensity = ca_data;
            cluster_traces(e,:) = ca_data;
            clear ca_data
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
 
        end
        
        activity_clusters.(['cluster_',num2str(k)]).mean_trace = mean(cluster_traces,1);
        xc = sum(x) / neighB_l;
        yc = sum(y) / neighB_l;
        activity_clusters.(['cluster_',num2str(k)]).barycenter(1,1:2) = [xc yc];
        
    end 
end
        
    %% 10 - Correlation matrix
    % We now build a correlation matrix based on mean fluorescence traces
    % of each cluster (calculated above)

    C = zeros(L,L);

    for k = 1:L
        for j = 1:L
            %Meantraces1 and 2
            Mt1 = activity_clusters.(['cluster_',num2str(k)]).mean_trace;
            Mt2 = activity_clusters.(['cluster_',num2str(j)]).mean_trace;
            X = corrcoef(Mt1,Mt2); 
            C(k,j) = X(1,2); %symmetrical: also equals to X(2,1)
        end
    end

    %% 11 - plots
    % A scatter plot with the registration template and a rasterplot 
    % are the two available options
    
    if logical(cl.spatialplot) 
        scatter_clusters(activity_clusters,frame_path) 
        
        figure; hold on
        imagesc(C), colorbar
        hold off 
    end
    
    if logical(cl.rasterplot)
        raster_clusters(activity_clusters); 
    end







end


function [cl] = check_empty(cl,field,default)

    if ~isfield(cl,field)
        cl.(field) = default;
    end


end
