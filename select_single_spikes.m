function [locs] = select_single_spikes(locs)

 global yes_events
 locs(find(yes_events == 0)) = [];
 
 