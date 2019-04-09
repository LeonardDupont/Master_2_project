function [] = manually_count_spikes(cnintensity)


    cnintensity = zero_and_max(cnintensity.'); 
    [pks,locs] = findpeaks(cnintensity,'MinPeakProminence',0.05);

     fig=figure('Name','Intensity = f(t)');  
     set(fig, 'KeyPressFcn',@keypress) %coupling the keypress with figure 
     global manual_spikes
     manual_spikes = [];
     peak = 1;
 

    function redraw() %will be executed each time until epoch == z
                clf
                hold on
                
                if peak < length(pks)
                    plot(cnintensity,'color','k')
                    plot(locs(peak),pks(peak) + 0.2,'*','color','blue')
                    axis tight
                
                
                    title(['Peak n� ', num2str(peak)])
                    xlabel('Points')
                    ylabel('Fluorescence')
                else
                    disp('Done.')
                end
                
     end
 
redraw()

       function keypress(~,evnt) %keypress reaction
            switch lower(evnt.Key)
                case 'rightarrow'
                    if peak<=length(pks)
                        manual_spikes(end+1) = locs(peak); %we keep its amplitude
                        peak = peak + 1;
                    end
                case 'leftarrow'
                    if peak<=length(pks)
                        peak = peak + 1;
                    end
                otherwise
                    return  %breaks it 
            end
            redraw()
       end

end