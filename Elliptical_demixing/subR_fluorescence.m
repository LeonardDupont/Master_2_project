function [neuropile] = subR_fluorescence(neuropile,inputfile)
% March 2019 - Carey lab - leonard.dupont@ens.fr
% .........................................................................
% This function takes the neuropile struct output from building_ellipses.m
% and builds fluorescence traces for the FISSA subregions using the calcium
% data specified in inputfile.
% .........................................................................
%
%  ----- INPUT ----------------------
%
%    neuropile      the output of building_ellipses, struct with fields
%
%    inputfile      either path or matrix
%  
%  ----- OUTPUT ---------------------
%
%    neuropile      struct with added field 'intensity'  
% .........................................................................


tic
disp('Reading tiff file...')
if isnumeric(inputfile)
    tiffdata = inputfile;
elseif ischar(inputfile)
    tiffdata = imread_tifflib(inputfile);
else
    error('Please attach compatible calcium data (path or matrix).') 
end
disp('Done.'), toc

[~,~,t] = size(tiffdata);
N = neuropile.n_cells;
nseg = neuropile.nseg;

disp('Now extracting fluorescence traces for neuropile subregions.')

for roi = 1:N
    np_mask_seg = neuropile.np_mask_seg{1,roi};
    
    if rem(roi,25) == 0 || roi == N || roi == 1
        disp(['roi ',num2str(roi),' out of ',num2str(N),'.']), toc
    end
        
    
    for seg = 1:nseg
        sub_seg = (np_mask_seg == seg); 
        sub_seg = uint16(sub_seg);
        seg_trace = zeros(1,t);
        
        for frame = 1:t
            a = tiffdata(:,:,frame) .* sub_seg;
            seg_trace(1,frame) = sum(a(:)) / length(find(a ~= 0)); 
        end    
        
        neuropile.intensity{seg,roi} = seg_trace; 
    end 
    
    center_roi = (np_mask_seg == nseg+1);
    center_roi = uint16(center_roi);
    center_trace = zeros(1,t);
    for frame = 1:t
        a = tiffdata(:,:,frame) .* center_roi;
        center_trace(1,frame) = sum(a(:)) / length(find(a~=0)); 
    end
    neuropile.intensity{seg+1,roi} = center_trace;
end

disp('Finished!'), toc

end