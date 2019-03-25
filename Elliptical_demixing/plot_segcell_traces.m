function [] = plot_segcell_traces(neuropile,cells)
%% March 2019 - Carey lab - leonard.dupont@ens.fr
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
% This function is to be used after applying subR_fluorescence.m in order
% to visualise both the mask of the neuropile and the calcium traces of the
% neuropile segments.
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
%
%  INPUT
%
%  neuropile     struct with fields, output of building_ellipses.m
%
%  cells         list of cell numbers to visualise (list or singulet)
%
%  OUTPUT
%
%  plot with normalised calcium data
%
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    
nseg = neuropile.nseg;

if neuropile.bkg
    nseg = nseg+1; 
    %if the background was added, there's an additional segment
end


N = length(cells);
cmap = parula(nseg+1); 

frames = length(neuropile.intensity{1,1}); 
t = linspace(1,frames/30,frames); % fs = 30

for i = 1:N
    cell = cells(i);
    % .   .   .   .   .    PLOTS THE MASK   .    .    .    .    .    .    . 
    figure, hold on
    subplot((nseg+1)+4,2,1:8) %first 8 slots are for the mask...
    imagesc(neuropile.np_mask_seg{1,cell});
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    start = 9; %...so we start from 9
    stop = 10;
    % .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 
    for seg = 1:nseg+1
        axh(seg) = subplot((nseg+1)+4,2,start:stop);
        ca = neuropile.intensity{seg,cell};
        ca = zero_and_max(ca); %normalise (linkedaxes = necessary)
        plot(t,ca,'color',cmap(seg,:))
        axis tight
        box off 
        set(gca,'ytick',[])
        set(gca,'TickLength',[0 0])
        
        start = stop + 1;
        stop = stop + 2;

        %now it's just details : ticks, titles and stuff
        if seg ~= nseg+1
            set(gca,'xtick',[])
        else
            ylabel('$\Delta F / F$')
            xlabel('time (s)')
        end  
        if neuropile.bkg
            if (seg == nseg) 
                title('soma')
            elseif (seg == nseg+1)
                title('background')
            else
                title(['seg ',num2str(seg)])
            end
        else
            if (seg == nseg+1)
                title('soma')
            else
                title(['seg ',num2str(seg)])
            end
        end
       
    end
    %finalising the visuals and linking the axis handles
    suptitle(['Cell ',num2str(cell)])
    linkaxes(axh,'xy')
    set(gcf,'Position',[500 500 400 1000])
end
