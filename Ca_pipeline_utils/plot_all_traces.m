function [] = plot_all_traces(cnintensity,spacing,fs)
% Simply plots all calcium recording traces of the detected ROIs on a
% single figure to get an idea of the diversity of signals

if nargin < 2
    spacing = 1000;
    fs = 30; %Hz
end

calcium_data = cnintensity; 
[x,y] = size(calcium_data);

clear back_to_zero
back_to_zero = zeros(x,y);

for roi=1:y
    back_to_zero(:,roi) = calcium_data(:,roi) - min(calcium_data(:,roi)); 
end

time = linspace(1,round(x/fs),x); 


figure, hold on
title('Traces of all detected ROIs in the FOV')
offset = 0;
positions = zeros(y,1); 
labels = cell(y,1);


for roi=1:y
    plot(time,back_to_zero(:,roi) + offset,'color','k')
    positions(roi,1) = offset;
    offset = offset + max(back_to_zero(:,roi)) + spacing; 
    labels(roi,1) = num2cell(roi);
end


x0 = 10; 
y0 = 10; 
width = 5500; 
height = 4000; 
set(gcf,'position',[x0,y0,width,height])
set(gca,'YTick',positions)
set(gca,'TickLength',[0.001,0])
set(gca,'yticklabel',labels)
box off 
axis tight
xlabel('Time (s)')
ylabel('ROIs')
hold off 