function [locs] = select_single_spikes(locs,yes_events)

 locs(find(yes_events == 0)) = [];
 
 