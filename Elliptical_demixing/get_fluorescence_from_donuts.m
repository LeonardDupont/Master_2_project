function [neuropile] = get_fluorescence_from_donuts(neuropile,inputfile,varargin)
% March 2019 - Carey lab - leonard.dupont@ens.fr
% ........................................................................
% This function uses the previously-defined, ellipse-shaped masks 
% surrounding the cells (building_ellipses.m) to extract a neuropile 
% fluorescence for each purkinje neuron. Each frame has a mean value.
% ........................................................................
%
%  --- INPUT -------
%
%   neuropile           struct with fields (see output of
%                       building_ellipses.m)
%   inputfile           path to the calcium recording (.tif)
%   varargin            weighed : 1 if using donut weights, 0 otherwise
%                       (default)
%
%  --- OUTPUT -------
%
%   neuropile          same as input, but added field 'intensity'
%
% .........................................................................

    ip = inputParser;
    ip.addParameter('weighed',0);
    parse(ip,varargin{:})
    weighed = logical(ip.Results.weighed);
 
%% Read tif

    tiffdata = imread_tifflib(inputfile);
    N = neuropile.n_cells;
    [h,w,t] = size(tiffdata);

%% Build mean neuropile fluorescence for each ROI (donut) and each frame 
    for roi = 1:N
        themask = uint16(neuropile.donutmask{1,roi});
        masked = uint16(zeros(h,w,t));
        
        for frame = 1:t
            a = tiffdata(:,:,frame).*themask; %convolute the mask and the FOV
            masked(:,:,frame) = a(:,:);
        end
        
        
        if weighed
            weights = neuropile.donutgradient{1,roi};
            weighed_mask = zeros(h,w,t);
            
            for frame = 1:t
                weighed_mask(:,:,frame) = masked(:,:,frame).* weights; %convolute the weights and
            end                               %the masked FOV
            
                    for frame = 1:t
                          theframe = weighed_mask(:,:,frame);
                          neuropile.intensity(frame,roi) = ...
                          sum(theframe(:))/length(find(theframe ~= 0)) ;
                    end  
            
        else
                    for frame = 1:t
                        theframe = masked(:,:,frame);
                        neuropile.intensity(frame,roi) = ...
                        sum(theframe(:))/length(find(theframe ~= 0)) ;
                    end
        end
    end
end