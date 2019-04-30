function [activity_clusters] = hybrid_clusteringDB(data,cl)
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
%        frame_path      path to registration template for spatialplot
%         rasterplot     logical, 1 if you want to plot (default : 0)
%          normalise     logical, 1 if you want to normalise (default : 0)
%        norm_method     either 'mean' or 'max' (default : [])
%               Z        minimum cluster size to be kept (discarded otherw)
%                        (default : 10)
%           gridbins     when using pixel-based clustering, number of bins
%                        to create in x and y directions (default : 15)
%            deconv0     set to 1 if you are using the 'deconvolution'
%                        structure post MLSpike (or other algorithm)
%
% -- OUTPUT -----------
%
%    activity_clusters      struct with fields
%                           
%                   
            

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
cl = check_empty(cl,'usedmdn',0); 
cl = check_empty(cl,'deconv0',0);
cl = check_empty(cl,'usedmdn',0);
cl = check_empty(cl,'Nmin',2);
cl = check_empty(cl,'Nmax',10);

if ~isstruct(data)
    disp('Using pixels as clustering units.')
    if ischar(data)
        name = data;
        clear data
        data = imread_tifflib(name);
        clear name
    end
    grdfov = grid_grouping(data,cl.gridbins);
    clear data
    N = (cl.gridbins)^2;
    [~,t] = size(grdfov.intensity);
    pixels = 1;
      
elseif isstruct(data)
    disp('Using rois as clustering units.')
    cn = data;
    N = cn.n_cells;
    pixels = 0;
    [~,t] = size(cn.intensity.');
    clear data
else
    error('Data and method do not match : check documentation. Aborting.')
end



%% Clustering units : regions of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~pixels 
    %% 1 - Normalisation 

    tic
    if cl.usedmdn && isfield(cn,'intensity_dm')
        traces = cn.intensity_dm.'; 
        if logical(cl.normalise)
            disp(' ------ 1 - Normalising the calcium data ------')
            switch cl.norm_method
                case 'mean'
                    zm_calcium_data = zero_and_mean(traces); 
                case 'max'
                    zm_calcium_data = zero_and_max(traces);
            end 
        else
            zm_calcium_data = traces; 
        end
    elseif cl.usedmdn && ~isfield(cn,'intensity_dm')
        warning('No field with demixed-denoised traces in cn. Using raw intensity instead')
        cl.usedmdn = 0;
    end
    if ~cl.usedmdn
        traces = cn.intensity.' ; 
        if logical(cl.normalise)
            disp(' ------ 1 - Normalising the calcium data ------')
            switch cl.norm_method
                case 'mean'
                    zm_calcium_data = zero_and_mean(traces); 
                case 'max'
                    zm_calcium_data = zero_and_max(traces);
            end 
        else
            zm_calcium_data = traces; 
        end
        
    end
    toc
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
    
    %% 6 - Hierarchical clustering 
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

    disp(' ------ 6 - Now merging using hierarchical clustering ------')

    distMAT = zeros(N,N);
    for k = 1:N
        R1 = neighbourhood.(['roi_',num2str(k)]).neighB;
        for j = 1:N
            R2 = neighbourhood.(['roi_',num2str(j)]).neighB;
            distMAT(k,j) = 1 - (length(intersect(R1,R2))/max(length(R1),length(R2)));
        end
    end
    
    y = squareform(distMAT);
    LK = linkage(y,'average');
    
    [Cvopt,Topt, ~] = ...
     find_cutoff_DB(LK,distMAT,'epsilon',1e-4,'Nmax',...
     cl.Nmax,'Nmin',cl.Nmin,'min_units',cl.Z);
    
    T = Topt;
    Nclust = max(T); 
    clear Topt
    
    disp(['Done! Number of clusters is ',num2str(Nclust),'.']); 
    toc
    
    %% 9 - Building output structure
    disp(' ------ 9 - Preparing output structure and correlation matrix -------')
%% Clustering units : rois %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ~pixels
    clear activity_clusters
    Treal = zeros(N,1);
    uu = 1;
    for k = 1:Nclust
   
        thecells = find(T==k);
        if length(thecells) > cl.Z
            Treal(thecells) = uu;
            activity_clusters.clusterregions{1,uu} = thecells;   
            l = length(thecells);
            cluster_traces = zeros(l,t);
            x = zeros(l,1);
            y = zeros(l,1);

            for j = 1:l
               roi = thecells(j);  
               % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
               posi = cn.centroid{1,roi};
               activity_clusters.centroid{uu,j} = posi;
               x(j) = posi(1,1);
               y(j) = posi(1,2);
               clear posi 
               % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
               ca_data = cn.intensity(:,roi);
               activity_clusters.intensity{uu,j} = ca_data;
               cluster_traces(j,:) = ca_data.';
               clear ca_data
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
               if logical(cl.deconv0)
                   spk = cn.spiketimes{1,roi};
                   activity_clusters.spiketimes{uu,j} = spk;
                   clear spk
               end
             end


            activity_clusters.mean_trace(:,uu) = mean(cluster_traces,1);
            xc = sum(x) / l;
            yc = sum(y) / l;
            activity_clusters.barycenter(1:2,uu) = [xc yc];
            uu = uu + 1;
        end
    end
    T = Treal; 
    Nclust = max(T);
    disp(['Number of real clusters is ',num2str(Nclust),'.'])
    
%% Clustering units : pixels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif pixels

  Treal = zeros(N,1);
  uu = 1;
  for k = 1:Nclust
      
     thecells = find(T==k);
     if length(thecells) > cl.Z
         Treal(thecells) = uu;
         activity_clusters.clusterregions{1,uu} = thecells;   
         l = length(thecells);
         cluster_traces = zeros(l,t);      
         x = zeros(l,1);
         y = zeros(l,1);

        for j = 1:l

            roi = thecells(j);
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            posi = grdfov.centroid{1,roi};
            activity_clusters.centroid{uu,j} = posi;
            x(e) = posi(1);
            y(e) = posi(2);
            clear posi
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            ca_data = grdfov.intensity(roi,:);
            activity_clusters.intensity{uu,j} = ca_data;
            cluster_traces(e,:) = ca_data;
            clear ca_data
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        end

        activity_clusters.mean_trace(:,uu) = mean(cluster_traces,1);
        xc = sum(x) / l;
        yc = sum(y) / l;
        activity_clusters.barycenter(1:2,uu) = [xc yc];
        uu = uu + 1;
     end
  end
    
   Nclust = max(T);
   T = Treal; 
   disp(['Number of real clusters is ',num2str(Nclust),'.'])
end


    clear reorganised
    reorganised = [];
    for k = 1:Nclust
        rois = find(T == k);
        l = length(rois);
        for i = 1:l
            reorganised(end+1) = rois(i);
        end
    end
    Nc = length(reorganised);
    orgdistMAT = zeros(Nc,Nc);
    for k = 1:Nc
        r1 = reorganised(k);
        for j = 1:Nc
            r2 = reorganised(j);
            orgdistMAT(k,j) = distMAT(r1,r2);
        end
    end

    figure, subplot(1,2,1), colormap('Gray'), imagesc(distMAT)
    subplot(1,2,2),  colormap('Gray'), imagesc(orgdistMAT)
    suptitle('Distance matrices before and after clustering'), hold off 


    %% 10 - Correlation matrix
    % We now build a correlation matrix based on mean fluorescence traces
    % of each cluster (calculated above)

    C = zeros(Nclust,Nclust);

    for k = 1:Nclust
        Mt1 = activity_clusters.mean_trace(:,k);
        for j = 1:Nclust
            %Meantraces1 and 2
            Mt2 = activity_clusters.mean_trace(:,j);
            X = corrcoef(Mt1,Mt2); 
            C(k,j) = X(1,2); %symmetrical: also equals to X(2,1)
        end
    end

    activity_clusters.corrMAT = C; 
    activity_clusters.distMAT = distMAT;
    activity_clusters.orgdistMAT = orgdistMAT;
    activity_clusters.Nclust = Nclust;
    activity_clusters.Ncells = length(find(T ~= 0));
    activity_clusters.treeCV = Cvopt; 
    
    %% 11 - plots
    % A scatter plot with the registration template and a rasterplot 
    % are the two available options
    
    if cl.spatialplot
        if isempty(cl.frame_path)
            warning('No FOV to scatter points on, please input path here.')
            cl.frame_path = input(prompt);
        end 
        
       if isstruct(cn)
            figure; hold on
            imh = imshow(cl.frame_path); hold on,
            title('Spatial distribution of clustered regions')
            cmap = jet(Nclust);
            S = size(cn.mask{1,1});
            for k = 1:Nclust
               regions = activity_clusters.clusterregions{1,k};
               L = length(regions);
               c = cmap(k,:);
               full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
               names = activity_clusters.clusterregions{1,k};
               for roi = 1:L
                 actual = regions(roi);
                 posi = activity_clusters.centroid{k,roi};
                 x = posi(1);
                 y = posi(2);
                 %text(x,y,num2str(names(roi)),'Color',c)
                 I = cn.mask{1,actual};
                 h = imshow(full); hold on 
                 set(h, 'AlphaData', I*0.35) , hold on
               end
            end
            % . . . . . . . . . . Lonely ROIS . . . . . . . . . . . . . . .
            k = 0;
            thecells = find(T == k);
            L = length(thecells);
            c = [1 1 1]; 
            for roi = 1:L
                 actual = thecells(roi);
                 posi = cn.centroid{1,actual};
                 x = posi(1);
                 y = posi(2);
                 %text(x,y,num2str(names(roi)),'Color',c)
                 I = cn.mask{1,actual};
                 h = imshow(full); hold on 
                 set(h, 'AlphaData', I*0.35) , hold on
            end
            hold off  
       else
           disp('Spatialplot not yet implemented for pixels in this function.')
       end
    
        figure; hold on
        imagesc(C), colorbar
        hold off 
    end
    
    if logical(cl.rasterplot)
        raster_clusters(activity_clusters); 
    end





disp('Clustering done!')

end


function [cl] = check_empty(cl,field,default)

    if ~isfield(cl,field)
        cl.(field) = default;
    end

end
