function [] = plot_spikes_w_speed(cnspikes,xc,speedlabels,tm_speedMs)

[~,Ncells] = size(cnspikes); 

%% 1 - Checking if tm_speedMs and cnintensity times match (in seconds)
% if huge discrepancy, then there is probably an issue and comparison is
% useless. this has been a recurrent problem because of timestamp gestion.

%Default
FsTM = 1/3.3333e-04; %Hz
FsIM = 30; %Hz
tTM = length(speedlabels)/FsTM;
tIM = length(cnspikes)/FsIM;
epsilon = 5; %seconds

if abs(tTM-tIM) > epsilon
    warning(['Oops! Treadmill time (',num2str(tTM),'s) and imaging time (',num2str(tIM),...
        's) do not seem to match... Di you use the right trial concatenation? Plotting nonetheless.'])
end

%% 2 - Rasterplot

height = 1;
spacing = 0.5;


for roi = 1:Ncells
     estimated = cnspikes(:,roi);
     spktimes = find(estimated==1)*FsTM/FsIM; %each round, we get the spike times 
     y = [(roi-1)*height + roi*spacing , roi*height + roi*spacing]; %we prepare 2 coordinates in y (we're gonna draw a vertical line: spike)

            for spk = 1:length(spktimes) %for each spike from this roi, we draw vertical lines at the right x coordinates
                x = [spktimes(spk), spktimes(spk)];
                plot(x,y,'color','k'), hold on
            end
end
   box off   

offset = roi*height + roi*spacing; 

%% 3 - Now plotting speed

offset = offset + 1000; 
position = offset;
timeTM = linspace(1,length(speedlabels),length(tm_speedMs));
plot(timeTM,tm_speedMs*300 + offset,'color','k','LineWidth',1)
offset = offset + max(tm_speedMs*300);

% This is just to label min speed, max speed and the speed trace
position = [position, (position+offset)/2 ,offset];
labels = {min(tm_speedMs),'speed (m.s^{-1})',max(tm_speedMs)};


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

%% 4 - Speed-based color filling of time bits

ncat = max(speedlabels);
colormap = hsv(ncat+1);


%  . . . . . . . . . Now we fill areas based on boundaries  . . . . . . . .
for k = 1:length(xc)-1
    xplus = xc(k+1);
    xminus = xc(k);
    label = speedlabels(xplus-1);
    color = colormap(label,:);
    
    xl = [xminus, xminus, xplus, xplus];
    yl = [0 offset offset 0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.15;
    a.EdgeColor = 'none';
end
hold on
xposi = linspace(1,length(speedlabels),10);
xlab = cell(10,1);
for k = 1:10
    xlab(k) = num2cell(round(xposi(k)/FsTM));
end

%  . . . . . . .  This to get a time axis in seconds  . . . . . . . . . . .
xticks(xposi)
xticklabels(xlab)

hold off 


