function [speedlabels] = assign_speedlabels(tm_speedMs,speedcats)
%% March 2019 - CareyLab - leonard.dupont@ens.fr
%  This function takes in speed as a function of time (tm_speedMs) and the
%  previously defined speedcategories (boundaries stored in speedcats) to
%  then label each datapoint with a category

    speedlabels = zeros(length(tm_speedMs),1);
    for j = 1:length(tm_speedMs)
        cat = 1;
        while tm_speedMs(j) > speedcats(cat)
            cat = cat + 1;
        end
        speedlabels(j) = cat;
    end

end