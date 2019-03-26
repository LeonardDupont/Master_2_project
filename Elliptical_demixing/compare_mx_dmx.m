function [] = compare_mx_dmx(cn,cells)
%% March 2019 - Carey lab - leonard.dupont@ens.fr
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
% This function is to be used after applying FISSA to visualise the
% demixing process in the rois specified in the 'cells' list. 
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
%
%  INPUT
%
%  cn            struct with fields, including intensity and intensity_dm
%
%  cells         list of cell numbers to visualise (list or singulet)
%
%  OUTPUT
%
%  plot with normalised calcium data
%
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

N = length(cells);
t = length(cn.intensity(:,1));
tt = linspace(0,t/30,t); %fs = 30Hz


for roi = 1:N
        cell = cells(roi);
        figure, hold on
        axh1 = subplot(2,3,1:3);
        mixed = zero_and_max(cn.intensity(:,cell).');
        plot(tt,mixed,'k')
        axis tight, box off
        title('mixed')
        axh2 = subplot(2,3,4:6); 
        demixed = zero_and_max(cn.intensity_dm(:,cell).');
        plot(tt,demixed,'k')
        axis tight, box off 
        title('demixed')
        xlabel('Time (s)')
        ylabel('\Delta F / F')
        linkaxes([axh1,axh2],'xy')
        suptitle(['Cell ',num2str(cell)])
        hold off
end