function [Dmin,neighbourhood,merge_events] = merge_closest(M,neighbourhood,N,merge_events,Dth)
%% February 2019 - Carey lab - leonard.dupont@ens.fr
%..........................................................................
% This function takes the similarity matrix M and merges close 
% neighbourhoods in the hybrid clustering procedure. 
% ........................................................................
% 
% -- INPUT ------------
%
%   M           matrix of distances between neighbourhoods as calculated by
%               distance_star.m, such that M(i,j) = distance between neighB
%               i and j.
%
%   nghbhood.   struct with fields as explained in hybrid clustering 
%
%   N           number of rois 
%
%   merge_evts  number of merge events that already occured since t=0
%
%   Dth         minimum distance (threshold) for merging
%
%
% -- OUTPUT ------------    
%
%   Dmin        minimal distance found in M
%
%   nghbhood.   updated structure after merging closest neighBs
%   
%   merge_evts  same as before, but +1 if Dmin was < Dth. 
% ........................................................................

  names = fieldnames(neighbourhood); %iscell == 1

  unwrapped = M(:);
  [Dmin,index] = min(unwrapped);
  
  if Dmin < Dth
  
      [i,j] = ind2sub(size(M),index);

       R1 = neighbourhood.(names{i}).neighB;
       R2 = neighbourhood.(names{j}).neighB;
       
       w1 = neighbourhood.(names{i}).weights;
       w2 = neighbourhood.(names{j}).weights;

       R3 = union(R1,R2);
       w3 = zeros(N,1);
       for k = 1:N
           w3(k) = (w1(k) + w2(k))/ 2; %new weight vector penalising unions
       end

       disp(['Now merging ',names{i},' and ',names{j},'.'])
       neighbourhood = rmfield(neighbourhood,names{i});
       neighbourhood = rmfield(neighbourhood,names{j});
       
       newname = ['roi_',num2str(N + merge_events)];
       neighbourhood.(newname).neighB = R3;
       neighbourhood.(newname).weights = w3;

       merge_events = merge_events + 1;
       
  else
      disp(['Minimum distance is bigger than threshold distance for merging. Number of clusters is ',num2str(length(names))])
      Dmin  = 987654321; 
  end
end