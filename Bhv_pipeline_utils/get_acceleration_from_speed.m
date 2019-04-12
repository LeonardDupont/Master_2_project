function [acceleration] = get_acceleration_from_speed(v,time)
% Takes speed and time, gives back the pointwise derivative (acceleration). 
% April 2019. 


L = length(v);
acceleration = zeros(L-1,1);
for t = 1:L-1
    acceleration(t) = ( v(t+1)-v(t) )/( time(t+1)-time(t) );
end
end