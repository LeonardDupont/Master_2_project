function [] = compare_spk2trace(cn,rois)


N = length(rois);
offset = 1;

for k = 1:N
    figure, hold on
    set(gcf, 'Position',  [100, 500, 1000, 200])
    trace = cn.intensity(:,rois(k)).';
    trace = zero_and_max(trace);
    ylim([0.4,max(trace + offset)+0.1])
    plot(trace + offset,'color','k'), hold on
    
    spikes = cn.spikes(:,rois(k)).'; 
    spktimes = find(spikes==1); 
    Nspk = length(spktimes);
    for i = 1:Nspk
        x = [spktimes(i),spktimes(i)];
        y = [0.5,0.8];
        plot(x,y,'color','b'), hold on
    end
    suptitle(['Cell ',num2str(rois(k))]); 
    
end