function   [pks,locs] = define_single_events(cn,rois)

    N = length(rois);
    disp('Concatenating traces.')
    t = length(cn.intensity(:,1)); 
    
    calcium_data = zeros(N*t,1); 
    for k = 1:N
        calcium_data(((k-1)*t+1):k*t,1) = cn.intensity(:,rois(k))-min(cn.intensity(:,rois(k))); 
    end
    
    %tt = linspace(1,N*t/30,N*t);
%% initial selection based on threshold

    
        [pks, locs,~,p] = findpeaks(calcium_data,'MinPeakProminence',0.05*max(calcium_data)); %theoretical minimum deltaF/F for one spike

        too_high = find(pks > 0.2*max(calcium_data));
        pks(too_high) = [];
        locs(too_high) = [];

        if isempty(pks)
            error('No single events detected in this ROI. Increase threshold or choose another one')
        end
 
 %% user-based selection
 
     disp('Initial selection of peaks performed successfully. Switching to user-based sorting. Press RIGHT ARROW to keep the peak (single-spike event), press LEFT ARROW otherwise.')

     fig=figure('Name','ROIs');  %overall figure name
     set(fig, 'KeyPressFcn',@keypress) %coupling the keypress with figure 
     global yes_events
     yes_events = [];
     peak = 1;
 

    function redraw() %will be executed each time until epoch == z
                clf
                hold on
                
                if peak < length(pks)
                    plot(calcium_data,'color','k')
                    plot(locs(peak),pks(peak) + 0.2,'*','color','blue')
                    axis tight
                
                
                    title(['Peak n° ', num2str(peak)])
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
                        yes_events(end+1) = p(peak); %we keep its amplitude
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
    