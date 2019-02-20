%% FEBRUARY 2019 - Carey lab
% This whole script is dedicated to an attempt of spatially clustering the
% cells based on their calcium activity. This all relies on mukamel rois,
% so there is no demixing yet. We will then try the pixel-based method,
% which should be much faster. 

%% normalise and bring back to baseline (0)

zm_calcium_data = zero_and_mean(calcium_data); 
[datapoints,rois] = size(zm_calcium_data);

%% Now we move onto clustering 
clear clusterresults
N_clusters = 5;
runs = 2500;


for i = 1:runs 
    idx = kmeans(zm_calcium_data.', N_clusters,'Distance','correlation'); %running k-means with
    for N = 1:N_clusters
        clusterresults.(['run_',num2str(i)]).(['c',num2str(N)]) = find(idx==N); 
    end
end 

%%

clear frequencies
ch_run = randsample(2500,1); %we randomly choose a run amongst all 


m = 6; %in each cluster of 'ch_run', we will pool 4 roi references
runs = numel(fieldnames(clusterresults)); %this is the number of k-means runs

for N = 1:N_clusters %now we will loop through each cluster
    regions = clusterresults.(['run_',num2str(ch_run)]).(['c',num2str(N)]);  %these are the rois from cluster N 
    l = length(regions); 
    samples = zeros(N_clusters,m);
    
    for sample = 1:m
        whichone = regions(randsample(l,1)); %we randomly choose a reference ROI
        samples(N,sample) = whichone; %and keep it there
        frequencies.(['clst_',num2str(N)]).(['sample_',num2str(sample)]).grouped = clusterresults.(['run_',num2str(ch_run)]).(['c',num2str(N)]);
    end


    %we now have m samples in cluster N. We need to find in which cluster they are in
    %the other (runs - 1) k-means... This could take some time. 
    run_nbs = linspace(1,runs,runs);
    run_nbs(ch_run) = []; 

    for sample = 1:m 
        for k = 1:(runs-1) %we used one for pooling

            run = run_nbs(k); %using the updated list (see above)
            c = 1;
            list = clusterresults.(['run_',num2str(run)]).(['c',num2str(c)]); 

            while isempty(find(list == samples(N,sample), 1)) %if the sample is not in this cluster...
                c = c + 1; %we try in the next one! 
                list = clusterresults.(['run_',num2str(run)]).(['c',num2str(c)]); 
            end

            konkat = cat(1,frequencies.(['clst_',num2str(N)]).(['sample_',num2str(sample)]).grouped, list);
            frequencies.(['clst_',num2str(N)]).(['sample_',num2str(sample)]).grouped = konkat;
            % and then we add the cluster from the tested run to our struct

        end
    end


    %for each of the m samples, we now have a huge vector that summarises
    %all of their local environments throughout the runs. Time to see if
    %these clusters are redundant or not. 
    [datapoints,rois] = size(zm_calcium_data);
    refs = linspace(1,rois,rois);

    for sample = 1:m

        stat_cluster = [];
        grouped = frequencies.(['clst_',num2str(N)]).(['sample_',num2str(sample)]).grouped;

        for roi = 1:rois

            probability = length(find(grouped == roi)) / runs; %number of appearances / number of runs = redundancy 

            if probability > 0.8
                stat_cluster(end+1) = roi; %if the roi is redundant enough our reference samples(m), we keep it! 
            end

        end

        if length(stat_cluster) > 10
            frequencies.(['clst_',num2str(N)]).(['sample_',num2str(sample)]).metacluster = stat_cluster; %this is the metacluster (stable rois)
        end
    end

    %at this point, we have metaclusters surrouding
    %each roi value (randomly chosen - m values) coming from cluster N: we estimated the stable
    %environment of each of these m rois (probability > 0.8). We now need
    %to fuse all of ours metaclusters IF THEY ARE REDUNDANT. We could do
    %that either by:
    %       (i) sheerly using the roi references in the metaclusters
    %       (ii) using the correlation between the mean calcium trace of
    %       each metacluster (Dombeck et al. (2009)). 
    % first things first : if there exists samples(mi) in one of the other
    % mj (i~=j) metaclusters, then the union is a metacluster. 

    frequencies.(['clst_',num2str(N)]).merged = [];
    for sample1 = 1:m
        for sample2 = 1:m
            
            if sample1 ~= sample2
                if isfield(frequencies.(['clst_',num2str(N)]).(['sample_',num2str(sample1)]) , 'metacluster') && isfield(frequencies.(['clst_',num2str(N)]).(['sample_',num2str(sample2)]) , 'metacluster')
                    meta1 = frequencies.(['clst_',num2str(N)]).(['sample_',num2str(sample1)]).metacluster;
                    meta2 = frequencies.(['clst_',num2str(N)]).(['sample_',num2str(sample2)]).metacluster;
                    isthere = intersect(meta1,meta2);

                    if ~isempty(isthere)
                        jointhem = union(meta1,meta2);
                        frequencies.(['clst_',num2str(N)]).merged = union(frequencies.(['clst_',num2str(N)]).merged,jointhem);
                    end
                else
                    warning(message(['No field "metacluster" for cluster ', num2str(N), ' and samples ', num2str(sample1), ' or ', num2str(sample2) '.'])
                end
 
            end

        end    
    end

    if numel(fieldnames(frequencies.(['clst_',num2str(N)]))) == 1
        warning(message(['No greater metacluster was created for cluster ',num2str(N), ' : intersection between sample metaclusters was null. Bad hint for cluster significance.']))
    end

end
%%
    %now we have merged a few metaclusters, we can check if some new
    %merging should happen. We'll use hierchical clustering at this point. 
    clear clusterstruct
   
    N = numel(fieldnames(frequencies));
    for cl = 1:N
         clusterstruct.(['cl_',num2str(cl)]) = frequencies.(['clst_',num2str(cl)]).merged ;
    end
    
    
    maximum = 50.9997;

    while maximum ~= 50.9997

        [maximum,max_c] = hierarchy(clusterstruct);
        names = fieldnames(clusterstruct);

        if maximum ~= 50.9997
               k = max_c(1);
               j = max_c(2);
               newmerge = intersect( clusterstruct.(names{k}) , clusterstruct.(names(j)) );
               clusterstruct = rmfield(clusterstruct,names{k});
               clusterstruct = rmfield(clusterstruct,names{j});
               clusterstruct.(['cl_',num2str(N+1)]) = newmerge;
        end
    end




function [maximum,max_c] = hierarchy(clusterstruct)

        maxi = 50.9997;
        max_c = [0 0];
        N = numel(fieldnames(clusterstruct));
        names = fieldnames(clusterstruct);
        intersection_M = zeros(N,N);
        for k = 1:N
            for j = 1:N
                group1 = clusterstruct.(names{k});
                group2 = clusterstruct.(names{j});

                intersection_M(k,j) = length(intersect(group1,group2)) / min(length(group1),length(group2)) ; 
            end
        end

        for k= 1:N
            for j = 1:N
                if intersection_M(k,j) > maxi
                   maxi = intersection_M(k,j); 
                   max_c = [k j];
                end 
            end
        end

        maximum = max(maxi,50.9997);

 end

%% get the calcium traces + the locations back



