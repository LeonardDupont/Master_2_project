function [fluolabels] = assign_fluolabels(fluodata,fluocats)
%% March 2019 - CareyLab - leonard.dupont@ens.fr
%  This function takes in speed as a function of time (tm_speedMs) and the
%  previously defined speedcategories (boundaries stored in speedcats) to
%  then label each datapoint with a category
    fluodata = zero_and_max(fluodata.'); 
    fluolabels = zeros(length(fluodata),1);
    for j = 1:length(fluodata)
        cat = 1;
        while fluodata(j) > fluocats(cat)
            cat = cat + 1;
        end
        fluolabels(j) = cat;
    end

end