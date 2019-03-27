function [] = plot_all_traces_w_speed(cnintensity,xc,speed_labels,ncat)
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

% ------------

colormap = jet(ncat+1);


% plot_dummies for legend
l1 = plot([NaN,NaN], 'color', colormap(1,:));
l2 = plot([NaN,NaN], 'color', colormap(2,:));
l3 = plot([NaN,NaN], 'color', colormap(3,:));
l4 = plot([NaN,NaN], 'color', colormap(4,:));
l5 = plot([NaN,NaN], 'color', colormap(5,:));
legend([l1, l2, l3, l4, l5], {'vlow','low','medium','fast','vfast'},'AutoUpdate','off');

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