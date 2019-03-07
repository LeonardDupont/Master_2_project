function [] = raster_clusters(activity_clusters)
% February 2019 - Carey lab - leonard.dupont@ens.fr
% .........................................................................
% This function uses the output of clustering_k_means and deconvolution 
% to create a colored rasterplot.
% .........................................................................
% 
% --- INPUT --------
%
%       activity_clusters     struct w\ fields, including spiketimes 
%
% --- OUTPUT --------
%
%       big, colored scatterplot
% .........................................................................

    N = numel(fieldnames(activity_clusters));
    spacing = 0.5;
    height = 1;
    
    middles = [];
    labels = {};
    
    rplot=figure('visible','off');
    hold on
    roi_counter = 0;
    
    for c = 1:N

        c_data = activity_clusters.(['cluster_',num2str(c)]);
        Nrois = numel(fieldnames(c_data.rois));
        roi_names = fieldnames(c_data.rois);
        colourmap = map(c,:);
        
        if c == 1
            check = c_data.rois(roi_names{1});
            if ~isfield(check,'spiketimes')
            warning('No spike times. Perform deconvolution.')
            error('Exiting rasterplot now.')
            end
        end

        for roi = 1:Nrois
            
            roi_counter = roi_counter + 1;
            spks = c_data.rois.(roi_names{roi}).spiketimes;
            number = strsplit(roi_names{roi},'_');
            number = number{2}; %just the number
            y = [(roi_counter-1)*height + roi_counter*spacing , roi_counter*height + roi_counter*spacing];
            middles(end+1) = roi_counter*height + roi_counter*spacing - 1/2*height; 
            labels{end+1} = number;
            
            for events = 1:length(spks)
                x = [spks(events),spks(events)];
                plot(x,y,'color',colourmap), hold on
            end
        end
    end
    
    set(groot,'CurrentFigure',rplot)
    xlabel('Time (s)')
    ylabel('Spatially-organised ROIs')
    set(gca,'TickLength',[0.001,0])
    set(gca,'Ytick',middles)
    set(gca,'YTickLabel',labels) %details to legend the rois over y 
    set(gca,'FontSize',30)
    box off
    axis tight
    hold off 
    set(rplot,'visible','on')


end