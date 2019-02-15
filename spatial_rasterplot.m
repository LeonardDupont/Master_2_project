function [spatial] = spatial_rasterplot(tiffimage,deconvolution,xbin)
% The purpose of this function is to spatially organise the
% previously-defined ROIs to do a first verification of the microzone
% hypothesis. We use the centroids. 

N = numel(fieldnames(deconvolution));
spacing = 0.5;
height = 1;

%% Show initial situation to user

h=figure('visible','off'); hold on
title('Spatial distribution of ROIs')

if ischar(tiffimage)
    frame = imread(tiffimage);  %read the path
elseif ismatrix(tiffimage)
    frame = tiffimage; %then we take the matrix
end

imshow(frame), hold on; %display fixed image of the calcium landscape

for roi = 1:N 
    x = deconvolution.(['roi_',num2str(roi)]).centroid(1,1);
    y = deconvolution.(['roi_',num2str(roi)]).centroid(1,2);
    scatter(x,y,36,[1 1 1],'filled') ; text(x,y,{num2str(roi)},'color','w')%plot the centroid position of each ROI with the ROI nb
end

set(h,'visible','on'); hold off 

%% Segregating segments in spatial bins

%creating spatial boundaries in terms of pixels
[~,w] = size(frame); %image width
pxMZ = floor(w/xbin); %pixelized MicroZones
aside = rem(w,xbin); 
boundaries = zeros(xbin,1);

for cat = 1:xbin
    boundaries(cat) = cat*pxMZ + (cat==xbin)*aside; %upper boundary
end 


%checking where each roi belongs to
for roi = 1:N
    x_pos = deconvolution.(['roi_',num2str(roi)]).centroid(1,1); %position in x
    region = 1;
    if x_pos < boundaries(region)
        region = 1;
    else
       while x_pos > boundaries(region)
        region = region + 1;
       end  
    end
    spatial.(['region_',num2str(region)]).(['roi_',num2str(roi)]).spiketimes = deconvolution.(['roi_',num2str(roi)]).spiketimes;
    spatial.(['region_',num2str(region)]).(['roi_',num2str(roi)]).centroid = deconvolution.(['roi_',num2str(roi)]).centroid;
end

%now we have a new struct called 'spatial'. Let's build the raster plot
%% Rasterplot with colormap + tiffimage with colormap

rplot=figure('visible','off');
hold on

rmap=figure('visible','off'); hold on
set(groot,'CurrentFigure',rmap)
imshow(frame), hold on;


cmap = colormap(parula(xbin));

roi_counter = 1; %we need to initialise a counter because of struct field names 

middles = zeros(N,1);
labels = cell(N,1);

for region = 1:xbin
    
    all_rois_reg = spatial.(['region_',num2str(region)]);
    subN = numel(fieldnames(all_rois_reg)); %number of ROIs belonging to region
    names = fieldnames(all_rois_reg); %their names (they don't go from 1 to subN: still initial numeration)
    
    for roi = 1:subN
        
        name = names{roi}; %the actual name (full)
        number = strsplit(name,'_');
        number = number{2}; %just the number
        x = all_rois_reg.(name).centroid(1,1);
        y = all_rois_reg.(name).centroid(1,2);
        set(groot,'CurrentFigure',rmap)
        scatter(x,y,36,cmap(region,:),'filled') ; text(x,y,number,'color',cmap(region,:)), hold on

        spktimes = all_rois_reg.(name).spiketimes;
        y = [(roi_counter-1)*height + roi_counter*spacing , roi_counter*height + roi_counter*spacing];
        middles(roi_counter) = roi_counter*height + roi_counter*spacing - 1/2*height; %get the middle coordinate of this line over y to label the ROI
        labels(roi_counter) = cellstr(number); %the tick label is the actual name, of course
        
        set(groot,'CurrentFigure',rplot)
        for spk = 1:length(spktimes) %for each spike from this roi, we draw vertical lines at the right x coordinates
                x = [spktimes(spk), spktimes(spk)];
                plot(x,y,'color',cmap(region,:)), hold on
        end
        
        roi_counter = roi_counter + 1;
    end
end

set(groot,'CurrentFigure',rplot)
xlabel('Time (s)')
ylabel('Spatially-organised ROIs')
set(gca,'TickLength',[0.001,0])
set(gca,'YTick',middles)
set(gca,'yticklabel',labels) %details to legend the rois over y 
box off 
axis tight
title(['Rasterplot of ', num2str(N), ' selected ROIs'])
hold off 

set(groot,'CurrentFigure',rmap)
imshow(frame), hold on;

set(rmap,'visible','on')
set(rplot,'visible','on')
    
end
    
    

    
    







