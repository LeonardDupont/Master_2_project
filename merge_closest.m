function [Dmin,neighbourhood,merge_events] = merge_closest(M,neighbourhood,rois,merge_events,Dth)
% 

  names = fieldnames(neighbourhood);

  unwrapped = M(:);
  [Dmin,index] = min(unwrapped);
  
  if Dmin < Dth
  
      [i,j] = ind2sub(size(M),index);

       R1 = neighbourhood.(names{i}).neighB;
       R2 = neighbourhood.(names{j}).neighB;
       
       w1 = neighbourhood.(names{i}).weights;
       w2 = neighbourhood.(names{j}).weights;

       R3 = union(R1,R2);
       w3 = zeros(rois,1);
       for k = 1:rois
           w3(k) = (w1(k) + w2(k))/ 2;
       end

       disp(['Now merging ',names{i},' and ',names{j},'.'])
       neighbourhood = rmfield(neighbourhood,names{i});
       neighbourhood = rmfield(neighbourhood,names{j});
       
       newname = ['roi_',num2str(rois + merge_events)];
       neighbourhood.(newname).neighB = R3;
       neighbourhood.(newname).weights = w3;

       merge_events = merge_events + 1;
       
  else
      disp(['Minimum distance is bigger than threshold distance for merging. Number of clusters is ',num2str(length(names))])
      Dmin  = 987654321; 
  end
end