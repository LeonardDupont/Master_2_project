function [] = show_trace_and_place(moviepath,calcium_data,roi_masks,list)
% This function is still slow, it might need to be optimised. 
% It plots a figure with the movie frames, roi masks and the corresponding
% fluorescence traces, all in different and corresponding colours. 
%
%  ------- INPUT ---------------
%
%   moviepath        path to the tiff stack to plot
%   calcium_data     the T * M matrix where T = frames (time) and M = rois
%   roi_masks        a cell with M entries, each containing a binary frame
%                    with 1s where the ROI is (see cn.mask)
%   list             vector of doubles corresponding to ROIs to be plot
%
%   ------ OUTPUT --------------
%
%   a figure with the movie playing as the fluorescence curves appear in
%   in the N = length(list) subplots under it. 
%
% February 2019 - Carey lab - leonard.dupont@ens.fr 



N = length(list); %number of rois to draw
cmap = parula(N); %preparing different colours with colormap

tiffmovie = bigread2(moviepath);
[~,~,Tf] = size(tiffmovie); 
Ts = linspace(1,Tf/30,Tf); %time in seconds if fs = 30Hz (default)

%%
used_data = zeros(Tf,N); %just for convenience, not essential
maxis = zeros(1,N);
minis = zeros(1,N);
for roi = 1:N
    used_data(:,roi) = calcium_data(:,list(roi));
    maxis(1,roi) = max(used_data(:,roi));
    minis(1,roi) = min(used_data(:,roi)); %useful to define ylims for each 
end                                       %subplot 

%%

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%turn it on if you want full screen

hAx(1)=subplot(6+3*N,1,1:6);
title('Field of view') %area for the movie frames

for roi = 1:N
        hAx(roi+1) = subplot(6 + 3*N,1,(6 + (roi-1)*3 + 1):(6 + roi*3));
        title(['roi_',num2str(list(roi))])
        ylabel('\Delta F / F_{0}')
        if roi == N
            xlabel(hAx(roi+1),'time (s)')
        end 
        ylim(hAx(roi+1),[minis(1,roi) maxis(1,roi)])
        box off 
end



for step = 1:Tf %and then we gradually plot the figure
    
    hAx(1) = subplot(6 + 3*N,1,1:6); % PLOT THE FRAME
    imshow(tiffmovie(:,:,step)), hold on
    title('Field of view')
    
    for roi = 1:N %PLOT BOUNDARY AGAIN ON FRAME 
       visboundaries(roi_masks{1,list(roi)},'color',cmap(roi,:),'Linewidth',0.05)
    end
    
    % now we update each subplot (fluorescence traces)
    
    for roi = 1:N
        hAx(roi+1) = subplot(6 + 3*N,1,(6 + (roi-1)*3 + 1):(6 + roi*3));
        plot(Ts(1:step),used_data(1:step,roi),'color',cmap(roi,:))
       
        if roi ~= N
            set(gca,'ytick',[])
            set(gca,'xtick',[])
        else %for the last one, we add labels and ticks
            set(gca,'ytick',[])
            set(gca,'TickLength',[0 0])
            ylabel(hAx(1+roi),'\Delta F / F_{0}')
            xlabel(hAx(1+roi),'time (s)')
        end
        xlim(hAx(roi+1),[1 Tf/30])
        ylim(hAx(roi+1),[minis(1,roi) maxis(1,roi)])
        box off
    end
    
    drawnow 
end