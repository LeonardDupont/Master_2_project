function [D] = distance_star(R1,R2,w1,w2)
% Calculates a weighed distance in the context of the clustering procedure.
% Please see the pdf file (clustering_analysis.pdf) or the script
% (clustering_k_means.m) for further explanations. 
%   
% INPUT
%        Ri     neighbourhood of region i (neighBi)
%        wi      weights for region i
%
%
% OUTPUT
%
%        D      distance between R1 and R2 as calculated hereunder
%
% Carey lab - February 2019

rois = length(w1);
card_inter = 0;
R1nR2 = intersect(R1,R2);

    for roi = 1:rois
        indicator = ~isempty( find(R1nR2 == roi) );
        card_inter = card_inter + w1(roi)*w2(roi)*indicator; 
    end

maxcard = max(length(R1),length(R2));

if card_inter ~= 0
    D = maxcard / card_inter; 
else
    D = 987654321;
end
    
end