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

    if isfield(activity_clusters,'lonely')
        K = numel(fieldnames(activity_clusters))-1;
        map = colormap(parula(K));

        h=figure('visible','off'); hold on
        ax = gca();
        hold(ax, 'on');
        title('Spatial distribution of clustered regions')

        for c = 1:K

            cluster = activity_clusters.(['cluster_',num2str(c)]).rois;
            Nrois = numel(fieldnames(cluster)); 
            roi_names = fieldnames(cluster);
            colourmap = map(c,:);

            for roi = 1:Nrois
                x = cluster.(roi_names{roi}).centroid(1);
                y = cluster.(roi_names{roi}).centroid(2);
                scatter(y,x,36,colourmap,'filled'), hold on
            end

            xb = activity_clusters.(['cluster_',num2str(c)]).barycenter(1,1);
            yb = activity_clusters.(['cluster_',num2str(c)]).barycenter(1,2);
            scatter(yb,xb,1,colourmap,'filled'), text(yb,xb,{num2str(c)},'color','w'); hold on

        end
        
        % . . . . . . LONELY REGIONS OF INTEREST . . . . . . . . . . . . . 
        cluster = activity_clusters.lonely.rois;
        Nrois = numel(fieldnames(cluster)); 
        roi_names = fieldnames(cluster);
        for roi = 1:Nrois
            x = cluster.(roi_names{roi}).centroid(1);
            y = cluster.(roi_names{roi}).centroid(2);
            scatter(y,x,36,'w','filled'), hold on
        end
        
        xb = activity_clusters.lonely.barycenter(1);
        yb = activity_clusters.lonely.barycenter(2);
        scatter(yb,xb,1,colourmap,'filled'), text(yb,xb,{num2str(c)},'color','w'); hold on
        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        
        imh = imshow(frame_path,'Parent',ax);
        uistack(imh,'bottom')
        hold(ax,'off')
        
        set(h,'visible','on'); hold off 
        
        
    else
        
        K = numel(fieldnames(activity_clusters));
        map = parula(K);

        h=figure('visible','off'); hold on
        ax = gca();
        hold(ax, 'on');
        title('Spatial distribution of clustered regions')

        for c = 1:K

            cluster = activity_clusters.(['cluster_',num2str(c)]).rois;
            Nrois = numel(fieldnames(cluster)); 
            roi_names = fieldnames(cluster);
            colourmap = map(c,:);

            for roi = 1:Nrois
                x = cluster.(roi_names{roi}).centroid(2);
                y = cluster.(roi_names{roi}).centroid(1);
                scatter(y,x,36,colourmap,'filled'), hold on
            end

            xb = activity_clusters.(['cluster_',num2str(c)]).barycenter(1);
            yb = activity_clusters.(['cluster_',num2str(c)]).barycenter(2);
            scatter(yb,xb,1,colourmap,'filled'), text(yb,xb,{num2str(c)},'color','w'); hold on

        end
        
        imh = imshow(frame_path,'Parent',ax);
        uistack(imh,'bottom')
        hold(ax,'off')
        
        set(h,'visible','on'); hold off 
        
        
    end

        

end