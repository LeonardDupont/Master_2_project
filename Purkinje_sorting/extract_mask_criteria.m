function [criteria] = extract_mask_criteria(cn)
%% March 2019 - Carey Lab - leonard.dupont@ens.fr
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
% This function extracts four features out of the cn masks in order to
% define if the segmented area is real or not. The four criteria are :
%
%     (i)  elongation       corresponds to the ratio of the 2 principal
%                           axes of the cell (mask). 
%                           'real' ones are elongated from what we know
%
%    (ii)  spread           the (%) of pixel surface in the FOV occupied
%                           by the cell. Too-small ones are usually removed
%    (iii) deviation        angle formed by the mask with the median field 
%                           angle (all rois detected) - based on the
%                           premise that most rois are real. 
%    (iv)  smoothness       regularity of the roi's surface, assessed by
%                           regionprops by defining a restricted ellipse
%                           based on the 2 first PCs and testing how much
%                           of the ROI is found inside of it. 
%
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
%
%   INPUT
%
%    cn        struct with fields (required : n_cells and mask)
%
%
%   OUTPUT
%
%    criteria  (n_cells * 4) double of observations (rows) and variables
%              (col)
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 


  N = cn.n_cells;
  Npx = cn.fov_height * cn.fov_width;
% -------------------------
  deviations = zeros(N,1);
  angles = zeros(N,1);
  spreads = zeros(N,1);
  elongations = zeros(N,1);
  smoothness = zeros(N,1);
% -------------------------
  
  for roi = 1:N 
      mask = cn.mask{1,roi}; 
      [xs,ys] = find(mask == 1);  
      points = [xs ys];

      COVM = cov(points);
      [V,D] = eigs(COVM);

      d = diag(D);
      e1 = [1,0];
      if sqrt(d(1)) > sqrt(d(2))
          elongations(roi) = sqrt(d(1)) / sqrt(d(2)); 
          v1 = V(1,1);
          v2 = V(1,2);
          angles(roi) = pi/2 + vector_angles2(e1,[v2,v1]);
      else
          elongations(roi) = sqrt(d(2)) / sqrt(d(1)); 
          v1 = V(2,1); %these are x = heights, so actually y in R^2
          v2 = V(2,2); %these are y = width, so actually x in R^2
          angles(roi) = pi/2 + vector_angles2(e1,[v2,v1]);
      end

      spreads(roi) = length(xs)*100/Npx; %in (%) 
      
      stats = regionprops(mask,'Solidity');
      smoothness(roi) = stats.Solidity; 
  end
  
  mdfield = median(angles);
  for roi = 1:N
      deviations(roi) = abs(mdfield - angles(roi));
  end
  
  criteria = [deviations,spreads,elongations,smoothness]; 
  
end
