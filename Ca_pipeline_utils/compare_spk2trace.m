function [] = compare_spk2trace(cn,rois,dm)


N = length(rois);
offset = 1;
if nargin < 3
    dm = 0;
end

if dm
    cnintensity = cn.intensity_dm;
else
    cnintensity = cn.intensity;
end
cmap = pink(4*N); 
for k = 1:N
    figure, hold on
    set(gcf, 'Position',  [100, 500, 1000, 200])
    trace = cnintensity(:,rois(k)).';
    trace = zero_and_max(trace);
    ylim([0.6,max(trace + offset)+0.1])
    xlim([0,length(trace)])
    plot(trace + offset,'color',cmap(k+N,:),'LineWidth',1.2), hold on
    
    spikes = cn.spikes(:,rois(k)).'; 
    spktimes = find(spikes==1); 
    Nspk = length(spktimes);
    for i = 1:Nspk
        x = [spktimes(i),spktimes(i)];
        y = [0.9,1];
        plot(x,y,'color',cmap(k+N,:),'LineWidth',1.2), hold on
    end
    title(['Cell ',num2str(rois(k))]); 
    
end