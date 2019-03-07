function [] = scatter_clusters(activity_clusters,frame_path)
% February 2019 - Carey lab - leonard.dupont@ens.fr
% .........................................................................
% This function uses the output of clustering_k_means to scatter cell
% clusters in space (using FOV frame). 
% .........................................................................
% 
% --- INPUT --------
%
%       activity_clusters     struct w\ fields  
%       frame_path            path to .tif file (FOV)
%
% --- OUTPUT --------
%
%       a plot with activity clusters spatially scattered in diffrt. colors
% .........................................................................


    N = numel(fieldnames(activity_clusters));
    map = colormap(parula(N));
    
    h=figure('visible','off'); hold on
    title('Spatial distribution of clustered regions')

    imshow(frame_path), hold on;

    for c = 1:N

        cluster = activity_clusters.(['cluster_',num2str(c)]);
        Nrois = numel(fieldnames(cluster.rois)) - 2; 
        %-2 : barycentre and mean trace
        roi_names = fieldnames(cluster.rois);
        colourmap = map(c,:);

        for roi = 1:Nrois
            x = cluster.rois.(roi_names{roi}).centroid(1,1);
            y = cluster.rois.(roi_names{roi}).centroid(1,2);
            scatter(y,x,36,colourmap,'filled'), hold on
        end
        
        x = cluster.barycenter(1,1);
        y = cluster.barycenter(1,2);
        scatter(y,x,1,colourmap,'filled'), text(x,y,{num2str(c)},'color','w'); hold on
        
    end

    set(h,'visible','on'); hold off 

end