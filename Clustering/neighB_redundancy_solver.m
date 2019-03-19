function [neighbourhood] = neighB_redundancy_solver(neighbourhood)
%% February 2019 - Carey lab - leonard.dupont@ens.fr
%..........................................................................
% This function finds elements (rois) that appear in more than one neighB
% and removes them from the neighB in which they have the smallest weight.
% ........................................................................
% 
% -- INPUT ------------
%
%  neighbourhood        struct with fields as described in
%                       hybrid_clustering
%
% -- OUTPUT ------------    
%
%   neighbourhood       struct with fields, but without redundancy btwn
%                       neighbourhoods
% .....................s...................................................


L = numel(fieldnames(neighbourhood));
names = fieldnames(neighbourhood);


    for k = 1:L
        for j = 1:L
             if k~=j
                 R1 = neighbourhood.(names{k}).neighB;
                 R2 = neighbourhood.(names{j}).neighB;

                 w1 = neighbourhood.(names{k}).weights;
                 w2 = neighbourhood.(names{j}).weights;

                 R1nR2 = intersect(R1,R2); %find common element (redundant)
                 inter = length(R1nR2); %how many of them are there?

                 if inter > 0
                     for i = 1:inter
                         redundant = R1nR2(i);
                         if w1(redundant) ~= w2(redundant)
                             closer = max(w1(redundant),w2(redundant));
                             
                             %if weight bigger in neighbourhood R1
                             if closer == w1(redundant)
                                 index = find(R2 == redundant);
                                 R2(index) = []; %remove it from R2
                                 w2(index) = 0;
                                 
                             %if weight bigger in R2    
                             elseif closer == w2(redundant)
                                 index = find(R1 == redundant);
                                 R1(index) = []; %remove it from R1
                                 w1(index) = 0;
                             end
                         end
                     end
                 end

                 neighbourhood.(names{k}).neighB = R1;
                 neighbourhood.(names{j}).neighB = R2;

                 neighbourhood.(names{k}).weights = w1;
                 neighbourhood.(names{j}).weights = w2;
             end
        end
    end



end