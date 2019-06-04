function [] = plot_all_traces_w_speed(cnintensity,xc,speed_labels,tm_speedMs,cnspikes)
% Plots all calcium traces specified in cnintensity on top of speed
% landscape with borders specified in xc and labels specified in
% speed_labels. If cnspikes is filled (optional argument), then the
% spiketrains are displayed below their corresponding traces (takes more
% time and GPU power). April 2019. 
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
%
%    INPUT
%
%    cnintensity    fluorescence intensity matrix 
%                   xdim : fluorescence, ydim : cells
%    xc             speedboundary limits (see get_speedregime_boundaries.m)
%    speed_labels   T * 1 vector of labels for the speed trace
%    tm_speedMs     speed trace, T * 1 vector of speed values (m.s-1)
%                   'treadmill_speedMeansmoothed'
%    cnspikes       spiketimes matrix
%                   xdim : spikes (0 & 1) , ydim : cells
%
%    OUTPUT
%
%    plot with color relating to  speed regimes
%
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
if nargin < 5
    wspikes = false;
else
    wspikes = true; %if cnspikes is an input
end
 

%% 1 - Checking if tm_speedMs and cnintensity times match (in seconds)
% if huge discrepancy, then there is probably an issue and comparison is
% useless. this has been a recurrent problem because of timestamp gestion.

%Default
FsTM = 1/3.3333e-04; %Hz
FsIM = 30; %Hz
tTM = length(speed_labels)/FsTM;
tIM = length(cnintensity)/FsIM;
epsilon = 5; %seconds

if abs(tTM-tIM) > epsilon
    warning(['Oops! Treadmill time (',num2str(tTM),'s) and imaging time (',num2str(tIM),...
        's) do not seem to match... Di you use the right trial concatenation? Plotting nonetheless.'])
end

%% 2 - Plotting fluorescence traces

timeIM = linspace(1,length(speed_labels),length(cnintensity));
[x,y] = size(cnintensity);
back_to_zero = zeros(x,y);
for roi=1:y
    back_to_zero(:,roi) = (cnintensity(:,roi) - min(cnintensity(:,roi)))/max(cnintensity(:,roi)); 
end



figure, hold on
imoffset = 0;
spkoffset = 0;
height = 0.05;
spacing = 0.01; 


for roi=1:y
    if wspikes
        spikes = cnspikes(:,roi);
        spktimes = find(spikes == 1) * FsTM/FsIM; 
        for k = 1:length(spktimes)
            xs = [spktimes(k), spktimes(k)];
            ys = [(roi-1)*height + roi*spacing + spkoffset , roi*height + roi*spacing + spkoffset];
            plot(xs,ys,'Color','k'), hold on
        end
        imoffset = imoffset + height; 
    end
    plot(timeIM,back_to_zero(:,roi) + imoffset,'color','k')
    imoffset = imoffset + max(back_to_zero(:,roi)) + wspikes*spacing; 
    spkoffset = spkoffset + max(back_to_zero(:,roi)); 
end

%% 3 - Now plotting speed

offset = imoffset + 0.2; 
position = offset;
timeTM = linspace(1,length(speed_labels),length(tm_speedMs));
plot(timeTM,tm_speedMs*max(back_to_zero(:))*20 + offset,'color','k','LineWidth',1)
offset = offset + max(tm_speedMs*max(back_to_zero(:))*20);

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
ylabel('\Delta F / F') %they occupy most of the space so OK 

%% 4 - Speed-based color filling of time bits

ncat = max(speed_labels);
clear colormap
colormap =  jet(ncat);  

% . . . . . . . . . . This is just to have a nice legend  . . . . . . . . . 
legend_spd = cell(ncat,1);
for i = 1:ncat
    l(i) = plot([NaN,NaN], 'color', colormap(i,:));
    spdlbl = ['speed_',num2str(i-1)];
    legend_spd(i) = cellstr(spdlbl); 
end
legend(l, legend_spd,'AutoUpdate','off');


%  . . . . . . . . . Now we fill areas based on boundaries  . . . . . . . .
for k = 1:length(xc)-1
    xplus = xc(k+1);
    xminus = xc(k);
    label = speed_labels(xplus-1);
    color = colormap(label,:);
    
    xl = [xminus, xminus, xplus, xplus];
    yl = [0 offset offset 0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.3;
    a.EdgeColor = 'none';
end
hold on
xposi = linspace(1,length(speed_labels),10);
xlab = cell(10,1);
for k = 1:10
    xlab(k) = num2cell(round(xposi(k)/FsTM));
end

%  . . . . . . .  This to get a time axis in seconds  . . . . . . . . . . .
xticks(xposi)
xticklabels(xlab)

hold off 



end