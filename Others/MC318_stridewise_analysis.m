%% MC318 - S4, stridewise analysis of calcium-imaging data
% In this script, we take the tracking data from S4 and perform a
% stride-wise analysis of the cell's responses to locomotion (proxy :
% calcium transients seen through GCaMP6f). 
%
% Several steps of post-processing : 
%
%       (1) Use the Mukamel masks which we have been taking in the previous
%           analysis and extract a session-wise fluorescence vector forall.
%
%       (2) Extract strides from the tracking, for each trial (in the right
%           order), build a cumulative session-wise time vector for
%           tracking and keep the cumulated time indices for start and stop
%           points of the considered step cycles using the tm.time vector. 
%           Strides will be sampled with 10 points. 
%
%       (3) Using the tdms time indices of strides from above and im.time,
%           find the corresponding start and stop fluorescence points. Be
%           over inclusive and take one more on each extremity. Then use
%           the maximal length time vector of them all to interpolate data
%           and align it to the right time (interp1d). 
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%   We will first do this for the mean fluorescence trace and then go for a
%   cell-resoluted analysis. We will plot data for the 4 paws. 



%% 1 - Extracting fluorescence from the whole session. 

wholeS4 = imread_tifflib(''); 
[h,w,t] = size(wholeS4); 
% . . . . . . replicating cns4 struct . . . . . . . . . . . . . . . . . . .
cnf.n_cells = cns4.n_cells;
cnf.fov_width = cns4.fov_width;
cnf.fov_height = cns4.fov_height;
cnf.mask = cns4.mask;
cnf.roi = cns4.roi;
cnf.roi_landscape = cns4.roi_landscape; 
cnf.centroid = cns4.centroid; 

% . . . . . . . . extracting sessionwise fluo . . . . . . . . . . . . . . .
for cc = 1:cnf.n_cells
    m = cnf.mask{1,cc}; 
    for fr = 1:t
        cnf.intensity(fr,:) = m .* wholeS4(:,:,fr); 
    end
end
cnf.mintensity = mean(cnf.intensity,2); %mean intensity over all rois 

% . . . . . . . . deconvolution  . . . . . . . . . . . . . . . . . . . . . 
ops.fs = 30;
ops.recomputeKernel =0;
ops.sensorTau = 0.7;
ops.estimateNeuropil = 0;
ops.deconvType = 'L0'; 
%ops.deconvType = 'OASIS';
threshold = 1000;
[~, ~, ~, cnfull] = get_spikes_from_calcium_traces(cnfull.intensity, ops, threshold, cnfull, []);


%% Extract strides

tracking_dir = ''; 

filePattern = fullfile(tracking_dir, '*.mat');
Ntrials = length(dir(filePattern));
trials = dir(filePattern); 



