function [] = plot_all_traces_w_speed(cnintensity,xc,speed_labels,tm_speedMs)
% Plots all calcium traces specified in cn_intensity on top of speed
% landscape with borders specified in xc and labels specified in
% speed_labels. 


calcium_data = cnintensity; 
[x,y] = size(calcium_data);
back_to_zero = zeros(x,y);
for roi=1:y
    back_to_zero(:,roi) = calcium_data(:,roi) - min(calcium_data(:,roi)); 
end
time = linspace(0,length(speed_labels),x);



figure, hold on
title('Traces of all detected ROIs in the FOV')
offset = 0;
positions = zeros(y,1); 
labels = cell(y,1);

for roi=1:y
    plot(time,back_to_zero(:,roi) + offset,'color','k')
    positions(roi,1) = offset;
    offset = offset + max(back_to_zero(:,roi)); 
    labels(roi,1) = num2cell(roi);
end

offset = offset + 1000; 
tt = linspace(1,length(speed_labels),length(tm_speedMs));
plot(tt,tm_speedMs*max(back_to_zero(:))*100 + offset,'color','k','LineWidth',1)
offset = offset + max(tm_speedMs*max(back_to_zero(:))*100);

x0 = 10; 
y0 = 10; 
width = 5500; 
height = 4000; 
set(gcf,'position',[x0,y0,width,height])
set(gca,'YTick',positions)
set(gca,'TickLength',[0,0])
%set(gca,'yticklabel',labels)
set(gca,'yticklabel',{})
box off 
axis tight
xlabel('Time (s)')
ylabel('ROIs')

% ------------
ncat = max(speed_labels);
colormap = hsv(ncat+1);

legend_spd = cell(ncat,1);
for i = 1:ncat
    l(i) = plot([NaN,NaN], 'color', colormap(i,:));
    spdlbl = ['speed_',num2str(i-1)];
    legend_spd(i) = cellstr(spdlbl); 
end
legend(l, legend_spd,'AutoUpdate','off');

for k = 1:length(xc)-1
    xplus = xc(k+1);
    xminus = xc(k);
    label = speed_labels(xplus-1);
    color = colormap(label,:);
    
    xl = [xminus, xminus, xplus, xplus];
    yl = [0 offset offset 0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.15;
    a.EdgeColor = 'none';
end


hold off 