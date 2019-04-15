function [] = fluo_x_accelerationpeaks(cnintensity,acc,colorgradient,xcacc)

FsTM = 1/3.333e-04;
timeIM = linspace(1,length(acc),length(cnintensity));
[x,y] = size(cnintensity);
back_to_zero = zeros(x,y);
for roi=1:y
    back_to_zero(:,roi) = cnintensity(:,roi) - min(cnintensity(:,roi)); 
end

figure, hold on
title('Fluorescence = f(t) with acceleration information')
offset = 0;
for roi=1:y
  plot(timeIM,back_to_zero(:,roi)+offset,'color','k'), hold on
  offset = offset + max(back_to_zero(:,roi)); 
end

offset = offset + 2e5; 
position = offset;
timeTM = linspace(1,length(acc),length(acc));
plot(timeTM,acc*max(back_to_zero(:))*100 + offset,'color','k','LineWidth',1), hold on
offset = offset + max(acc*max(back_to_zero(:))*100);


for k = 1:length(xcacc)-1
    xplus = xcacc(k+1);
    xminus = xcacc(k);
    color = colorgradient(xplus-1,:);
    
    xl = [xminus, xminus, xplus, xplus];
    yl = [0 offset offset 0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.2;
    a.EdgeColor = 'none';
end



% This is just to label min speed, max speed and the speed trace
position = (position+offset)/2;
labels = {'acceleration (m.s^{-2})'};

x0 = 10; 
y0 = 10; 
width = 5500; 
height = 4000; 
set(gcf,'position',[x0,y0,width,height])
set(gca,'YTick',position)
set(gca,'TickLength',[0,0])
set(gca,'yticklabel',labels)
box off 
axis tight
xlabel('Time (s)')
ylabel('ROI fluorescence') %they occupy most of the space so OK 

accmaps = [0 0.3 0.9 ; 0.9 0.9 0.9 ;1 0.4 0.3]; 
accnames = {'decceleration','constant speed','acceleration'};
for i = 1:3
    l(i) = plot([NaN,NaN], 'color', accmaps(i,:));
end
legend(l, accnames,'AutoUpdate','off');


hold on
xposi = linspace(1,length(acc),10);
xlab = cell(10,1);
for k = 1:10
    xlab(k) = num2cell(round(xposi(k)/FsTM));
end


%  . . . . . . .  This to get a time axis in seconds  . . . . . . . . . . .
xticks(xposi)
xticklabels(xlab)

hold off 

end