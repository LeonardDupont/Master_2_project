function [zm_calcium_data] = zero_and_mean(calcium_data)
% produces a deltaF over F where centering occurs using the trace minimum
% value and normalisation occurs using the mean.
%
% INPUT
%
%    calcium_data       xdim = observations, ydim = variables (cells)
%
% OUTPUT
%
%    zm_calcium_data    zero + mean (zm) calcium data 
%
% Carey lab - February 2019

[x,y] = size(calcium_data);
zm_calcium_data = zeros(x,y);

for roi = 1:y
    mini = min(calcium_data(:,roi));
    moyenne = mean(calcium_data(:,roi)); 
    zm_calcium_data(:,roi) = (calcium_data(:,roi) - mini) / moyenne; 
end

