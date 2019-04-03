%% CareyLab - 2019 - Working Script 2
% This script uses the roi data extracted in WS2 to make the fluorescence
% signals revelant and draw info. from them. 
%
%       (A)  Elliptical demixing of ROIs 
%       (B)  Calcium signal deconvolution 
%       (C)  Clustering algorithms 
%       (D)  ISI histograms, rasterplots
%       (PL) Plotting tools
%
% In order to easily reproduce the performed analysis, this script should
% never be used as such. It should be copied and pasted in a new .m file
% that will only be used for the analysis of a given session. In other
% words, each session analysis should have its own WS script folder
% accompanying the data (see Z).
%
% Contributors : Hugo Marques, Leonard Dupont.
% leonard.dupont@ens.fr
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

%% (Z) Analysis information

WSanalysis.dd =  date; %date of analysis, dday
WSanalysis.mouse = 'MC___'; %please complete
WSanalysis.session = 'SX';


%% (A) Non-negative matrix factorisation to demix 1P calcium signals
% All explanations are found within the used functions. 

% A1 - This if you want to include background in demixing :
% eroding might be an option if it's a large area  
msk = cn_bkg.mask{1,1};
msk = imerode(msk,strel('disk',20));
cn_bkg.mask{1,1} = msk; 


% A2 - Build ellipsoids around cell bodies - create be options 
nseg = 4; 
bkg = true;
bkgmask = cn_bkg.mask{1,1};
l = 7;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
neuropile = ...
    building_ellipses(cn,'graphics',0,'segment_ellipse',1,'nseg',nseg,'bkg',bkg,'bkgmask',bkgmask,'l',l);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% A3 - Extract fluorescence values corresponding to the elliptic masks
inputfile = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51.tif';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
neuropile = subR_fluorescence(neuropile,inputfile); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% A4 - Demixing using NNMF
cn = FISSAdemix_from_ellipses(neuropile,cn); 


%% (B) - Deconvolution of calcium signals
% According to Pachitariu and Kenneth (2018), best results are provided by
% nonconstrained positive kernels such as OASIS in the case of data with no
% prior knowledge. However in our case, it seems that the L0-norm constrain
% works well.


% B1 - defining single events to adjust the threshold
% Choose rois that seem to display single events, list them and define the
% spikes using the functions hereunder.
plot_all_traces(cn.intensity)
se_traces = [23 24]; 
define_single_events(cn,se_traces)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
eventsamp = getGlobal_yesevents; 
se_amp = mean(eventsamp); %mean single-event amplitude
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% B2 - deconvolution using L0 norm (what works best right now)
ops.fs = 30;
ops.recomputeKernel =0;
ops.sensorTau = 0.7;
ops.estimateNeuropil = 0;
ops.deconvType = 'L0'; 
%ops.deconvType = 'OASIS';
threshold = se_amp;
[~, ~, ~, cn] = get_spikes_from_calcium_traces(cn.intensity, ops, threshold, cn, []);
%[~, ~, ~, cn] = get_spikes_from_calcium_traces(cn.intensity_dm, ops, threshold, cn, []);


% B1prime - using Thomas Deneux's MLSpike algorithm (needs rework) 
edit WS2_annex_MLSpike.m 



%% (C) - Clustering algorithms 
% Two available options:
%       (C1) Hierarchical clustering based on synchrony (needs
%       deconvolution prior to use).
%       (C2) Improved k-means clustering based on fluorescence of either
%       rois or pixels. 


% C1 - Synchrony
clear opts
opts.Nmax = 10;
opts.Nmin = 2;
opts.epsilon = 10e-4;

clear grphcs
grphcs.dendrogram = 1;
grphcs.orgdistMAT = 1;
grphcs.orgchanceMAT = 1;
grphcs.porgchanceMAT = 1;
grphcs.clusteredlandscape = 1;

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
synchresults = synchrony_clustering(cn,opts,grphcs); 
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


% C2 - Fluorescence
clear cl
cl.normalise = 1;
cl.K = 4;
cl.Dth = 1/0.4; 
cl.p = 0.6; 
cl.spatialplot = 1;
cl.usedmdn = 0; 
cl.gridbins = 15; 
cl.frame_path = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/registration_template.tif';
video = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51.tif';

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
[activity_clusters,C] = hybrid_clustering(video,cl);
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .




%% (D) - Additional analysis 

% D1 - Build ISI histogram
% See function description, can be used on all cells or a selection, etc.
build_ISI_histo(cn)

% D2 - Rasterplot
rasterplot(cn.spikes.')


%% (PL) - Plotting functions

% PL1 - to visualise the ellipse, the fluorescence of the neuropile
% segments and compare it to the center ROI.
cells = [1,2,3]; 
plot_segcell_traces(neuropile,cells)

% PL2 - after NNMF, FISSA demixing : compare ROI fluorescence before and
% after the demixing protocol. 
cells = [1,2,3]; 
compare_mx_dmx(cn,cells)

% PL3 - compare the fluorescence trace and the estimated spiketrain after
% deconvolution. 
rois = [1 2 3];
compare_spk2trace(cn,rois)

% PL4 - more global
plot_all_traces(cn.intensity)
plot_all_traces(cn.intensity_dm)


%% Consider saving the workspace?

wkspacename = ['WS1_workspace',WSanalysis.dd,WSanalysis.mouse,WSanalysis.session]; 
save(wkspacename)




