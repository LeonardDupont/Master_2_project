function [xc] = get_speedregime_boundaries(speedlabels)

    xc = [];
    label = 0;
    for p = 1:length(speedlabels)
        nlabel = speedlabels(p);
        if nlabel ~= label
            xc(end+1) = p;
            label = nlabel;
        end
    end
    xc(end+1) = length(speedlabels);

end