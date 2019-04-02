function [speedcats] = define_speedcats(tm_speedMs,ncat)
%% March 2019 - CareyLab - leonard.dupont@ens.fr
% If tm_speedMs is the speed as a function of time, ncat the number of
% speed categories that the user wants to create, then this function
% outputs speedcats, a vector defining the speed regime boundaries.
%
% ex: the speed ranges from 0 to 1 and ncat = 4, then we will output
% speedcats = [ 0  , 0.25  , 0.5  , 1]


    tm_speedMs(tm_speedMs < 0) = 0; %bringing back to 0
    maxs = max(tm_speedMs);
    mins = min(tm_speedMs);


    diffcat = (maxs - mins)/(ncat) ;
    speedcats = zeros(ncat+1,1);
    for i = 1:ncat+1
        speedcats(i) = mins + (i-1) * diffcat;
    end

end