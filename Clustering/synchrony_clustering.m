function [results] = synchrony_clustering(cn,opts,grphcs)
%% April 2019 - CareyLab - leonard.dupont@ens.fr
% ......................................................................... 
% This function is a wrapper that uses the synchrony of action potentials
% to cluster cells based on common activity pattern in time. To do so, it
% uses a measure of synchrony (defined hereunder) and hierarchical
% clustering. To make sure that the grouping is statistically relevant, we
% compute -for each pair of cells- the probability to recapitulate their
% synchrony value solely by chance using a Poisson process. Besides, we
% optimise the cutoff value of our tree according to an adapted version of
% the Davies-Bouldin criterion.
% Advice from Jorge Ramirez (PhD) was provided.  
% .........................................................................
% 
%   INPUT
%
%   cn         struct with (necessary) fields :
%              {n_cells, spikes, centroid, mask}
%
%   opts       struct of options with fields:
%                .framepath   @char        path to registration template
%                .Nmax         @double     maximal number of clusters
%                .Nmin         @double     minimal number of clusters
%                .epsilon      @double     step used in find_cutoff_DB.m
%                .min_units    @double     minimal nb of units to count cl
%                .wsize        @double     window size (frames) to consider
%                                          synchrony
%                .computechance @logical   if you want to compute the
%                                          chanceMAT
%
%   grphcs     struct of options with fields (all boolean)
%                .dendrogram
%                .distMAT
%                .chanceMAT
%                .orgdistMAT
%                .orgchanceMAT
%                .porgchanceMAT (orgchanceMAT < 0.05 - logicals)
%                .clusteredlandscape
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%
%   OUTPUT
%
%   results    struct with fields:
%                .synchMAT     @N*Nmat      matrix with synchrony values 
%                .distMAT      @N*Nmat      1 - synchMAT : distance matrix
%                .chanceMAT    @N*Nmat      matrix with probabilities to 
%                                           recapitulate synchrony values
%                                           of synchMAT by chance
%                .orgchanceMAT              chanceMAT organised by clusters
%                .orgdistMAT                distMAT organised by clusters
%                .clusteredrois             struct with fields containing
%                                           clustering results:
%                                          {clusterregions,centroids,masks}
%                .Cvopt       @double      optimal (and used) cutoff value
%                .T           @N*1vect     clustering indexes for the cells
%
%   figures    depending on the graphics options
%
% .........................................................................

N = cn.n_cells;
opts = check_empty(opts,'Nmax',10);
opts = check_empty(opts,'Nmin',2);
opts = check_empty(opts,'epsilon',1e-4);
opts = check_empty(opts,'framepath',[]);
opts = check_empty(opts,'min_units',3); 
opts = check_empty(opts,'wsize',2); 
opts = check_empty(opts,'computechance',1); 

grphcs = check_empty(grphcs,'dendrogram',0);
grphcs = check_empty(grphcs,'distMAT',0);
grphcs = check_empty(grphcs,'chanceMAT',0);
grphcs = check_empty(grphcs,'orgdistMAT',0);
grphcs = check_empty(grphcs,'orgchanceMAT',0);
grphcs = check_empty(grphcs,'porgchanceMAT',0);
grphcs = check_empty(grphcs,'clusteredlandscape',0);

%% 1 - Estimating paired synchrony values based on spike trains
disp('1 - Building synchrony and chance matrices out of spiketimes.')

synchMAT = zeros(N,N);
chanceMAT = zeros(N,N);
tic
wsize = opts.wsize;
cnspk = cn.spikes; 
computechance = opts.computechance; 
for i = 1:N
    x = cnspk(:,i).';
    parfor j = 1:N
        y = cnspk(:,j).';
        [Q,~,~,taumin] = est_trace_synchrony(x,y,wsize);
        synchMAT(i,j) = Q;
        distMAT(i,j) = 1 - synchMAT(i,j);
        if computechance
            chanceMAT(i,j) = est_synchrony_reliability(x,y,taumin); 
        end
    end
    percent = i*100/N; 
    if rem(floor(percent),5) == 0
        disp([num2str(floor(percent)),'% done.'])
    end
end

disp('Done!')
toc
 

%% 2 - Building a tree with nodes from the distance matrix
% For this we use the linkage function of MathWorks clustering toolbox.
% - squareform unwraps distMAT into a compatible vector.
% - 'average' is the method to estimate cluster distances (UPGMA)
disp('2 - Building tree.')

y = squareform(distMAT);
Z = linkage(y,'average');

disp('Done!')
toc

%% 3 - Finding the optimal cutoff value using Davies-Bouldin and clustering
% All explanations can be found in the function below.
disp('3 - Determining optimal number of clusters using Davies-Bouldin.')

[Cvopt,Topt, ~] =...
     find_cutoff_DB(Z,distMAT,'epsilon',opts.epsilon,'Nmax',...
     opts.Nmax,'Nmin',opts.Nmin,'min_units',opts.min_units);
 
if grphcs.dendrogram
    H = dendrogram(Z,'Orientation','left','ColorThreshold','default');
    vline(Cvopt), hold off
end

T = Topt;
Nclust = max(T); 
clear Topt

disp(['Done! Number of clusters is ',num2str(Nclust),'.']); 
toc
%% 4 - Organising data, matrices and preparing output plots
disp('4 - Preparing plots and output structure.')


N = cn.n_cells; 
clear activity_cluster
clear reorganised
reorganised = [];
for k = 1:Nclust
    rois = find(T == k);
    l = length(rois);
    activity_cluster.clusterregions{1,k} = rois;
    for i = 1:l
        reorganised(end+1) = rois(i);
        activity_cluster.centroid{k,i} = cn.centroid{1,rois(i)};
        activity_cluster.mask{k,i} = cn.mask{1,rois(i)};
    end

end

orgdistMAT = zeros(N,N);
orgchanceMAT = zeros(N,N);
for k = 1:N
    r1 = reorganised(k);
    for j = 1:N
        r2 = reorganised(j);
        orgdistMAT(k,j) = distMAT(r1,r2);
        orgchanceMAT(k,j) = chanceMAT(r1,r2);
    end
end

if grphcs.orgdistMAT
    figure, imagesc(orgdistMAT), title('Distance matrix with clustering-based grouping'), colorbar
end
if grphcs.distMAT
    figure, imagesc(distMAT), title('Distance matrix with no grouping'), colorbar
end
if grphcs.chanceMAT
    figure, imagesc(chanceMAT), title('Matrix showing P_{chance} with no grouping'), colorbar
end
if grphcs.orgchanceMAT
    figure, imagesc(orgchanceMAT), title('Matrix showing P_{chance} with clustering-based grouping'), colorbar
end
if grphcs.porgchanceMAT
    orgchanceMAT(isnan(orgchanceMAT)) = 0;
    porgchanceMAT = orgchanceMAT < 0.05;
    figure, imshow(porgchanceMAT); 
end

if grphcs.clusteredlandscape
     imh = imshow(opts.framepath); hold on,
     title('Spatial distribution of clustered regions')
     cmap = jet(Nclust);
     S = size(cn.mask{1,1});
     for k = 1:Nclust
         L = length(activity_cluster.clusterregions{1,k});
         c = cmap(k,:);
         full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
         names = activity_cluster.clusterregions{1,k};
         for roi = 1:L
             posi = activity_cluster.centroid{k,roi};
             x = posi(1);
             y = posi(2);
             text(x,y,num2str(names(roi)),'Color',c)
             I = activity_cluster.mask{k,roi};
             h = imshow(full); hold on 
             set(h, 'AlphaData', I*0.35) , hold on
         end
     end
     hold off  
end


% . . . . . . . . . . . . . . . OUTPUT STRUCT . . . . . . . . . . . . . . . 
results.synchMAT = synchMAT;
results.distMAT = distMAT;
results.chanceMAT = chanceMAT;
results.orgdistMAT = orgdistMAT;
results.orgchanceMAT = orgchanceMAT;
results.clusteredrois = activity_cluster;
results.Cvopt = Cvopt;
results.T = T;
% . . . . . . . . . . . . . . . . . . . . . . . . .. . . . . . . . . . . . 

disp('Exiting. Clustering done!')


    function [strct] = check_empty(strct,field,default)

        if ~isfield(strct,field)
            strct.(field) = default;
        end

    end

end



