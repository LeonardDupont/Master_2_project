function [fluocats] = define_fluocats(fluodata,ncat)
%% March 2019 - CareyLab - leonard.dupont@ens.fr
% Same as define_speedcats.m, but with fluorescence. 


    fluodata = zero_and_max(fluodata);
    maxs = max(fluodata);
    mins = min(fluodata);


    diffcat = (maxs - mins)/(ncat) ; %mins is 0, fluorescence rarely is 0
    fluocats = zeros(ncat,1);
    for i = 1:ncat
        fluocats(i) = i* diffcat;
    end

end