function   [pks,locs] = define_single_events(calcium_data)

%% initial selection based on threshold

    if mean(calcium_data) ~= 1
        calcium_data = calcium_data/mean(calcium_data);
    end
    
    [pks, locs] = findpeaks(calcium_data,'MinPeakProminence',0.0341*max(calcium_data)); %theoretical minimum deltaF/F for one spike
    
    too_high = find(pks > 0.7*max(calcium_data));
    pks(too_high) = [];
    locs(too_high) = [];
    
    if isempty(pks)
        error('No single events detected in this ROI. Increase threshold or choose another one')
    end
 
 %% user-based selection
 
 disp('Initial selection of peaks performed successfully. Switching to user-based sorting. Press RIGHT ARROW to keep the peak (single-spike event), press LEFT ARROW otherwise.')
  
 fig=figure('Name','ROIs');  %overall figure name
 set(fig, 'KeyPressFcn',@keypress) %coupling the keypress with figure 
 yes_events = zeros(length(pks),1);
 peak = 1;
 

    function redraw() %will be executed each time until epoch == z
                clf
                hold on
                
                if peak < length(pks)
                    plot(calcium_data)
                    plot(locs(peak),pks(peak) + 0.2,'*','color','red')
                    axis tight
                
                
                    title(['Peak n° ', num2str(peak)])
                    xlabel('time (ms)')
                    ylabel('Normalised fluorescence')
                else
                    disp('Done.')
                end
                
     end
 
redraw()

       function keypress(~,evnt) %keypress reaction
            switch lower(evnt.Key)
                case 'rightarrow'
                    if peak<=length(pks)
                        yes_events(peak) = 1; %we keep it
                        peak = peak + 1;
                    end
                case 'leftarrow'
                    if peak<=length(pks)
                        yes_events(peak) = 0; %we do not
                        peak = peak + 1;
                    end
                otherwise
                    return  %breaks it 
            end
            redraw()
       end
   
global yes_events

end
    