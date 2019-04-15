function [] = rasterplot(data, varargin)
%This function is used in the context of the calcium imaging pipeline. 
% February 2019 - Carey lab (LD)
% INPUT
%
%    data       n * p binary matrix (timepoints w\ 1 if spike)
%               n = nb of cells, p = datapoints
%
% OUTPUT
% 
%     the function outputs a regular raster plot with ROIs over y and time in
%     seconds over x. 

ip = inputParser;
ip.addParameter('height',1);
ip.addParameter('spacing',0.5);
ip.addParameter('color','k');

parse(ip,varargin{:});

height = ip.Results.height;
spacing = ip.Results.spacing;
colour = ip.Results.color;

%%
    [N,~] = size(data); 
    
    figure, hold on
    xlabel('Time (s)')
    ylabel('ROI')
    middles = zeros(N,1);
    labels = cell(N,1);
    
    for roi = 1:N
        
        estimated = data(roi,:);
        spktimes = find(estimated==1); %each round, we get the spike times 
        y = [(roi-1)*height + roi*spacing , roi*height + roi*spacing]; %we prepare 2 coordinates in y (we're gonna draw a vertical line: spike)
        middles(roi) = roi*height + roi*spacing - 1/2*height; %get the middle coordinate of this line over y to label the ROI
        labels(roi) = num2cell(roi);
        
            for spk = 1:length(spktimes) %for each spike from this roi, we draw vertical lines at the right x coordinates
                x = [spktimes(spk), spktimes(spk)];
                plot(x,y,'color',colour)
            end
        
    end
    
    set(gca,'TickLength',[0.001,0])
    set(gca,'YTick',middles)
    set(gca,'yticklabel',labels) %details to legend the rois over y 
    box off 
    axis tight
    title(['Rasterplot of ', num2str(N), ' selected ROIs'])

    
end