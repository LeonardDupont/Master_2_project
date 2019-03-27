function [zm_calcium_data] = zero_and_max(calcium_data)
% produces a deltaF over F where centering occurs using the trace minimum
% value and normalisation occurs using the max.
%
% INPUT
%
%    calcium_data       xdim = cells, ydim = time
%
% OUTPUT
%
%    zm_calcium_data    zero + max (zm) calcium data 
%
% Carey lab - February 2019

[x,y] = size(calcium_data);
zm_calcium_data = zeros(x,y);

for roi = 1:x
    mini = min(calcium_data(roi,:));
    maxi = max(calcium_data(roi,:)-mini); 
    zm_calcium_data(roi,:) = (calcium_data(roi,:) - mini) / maxi; 
end

