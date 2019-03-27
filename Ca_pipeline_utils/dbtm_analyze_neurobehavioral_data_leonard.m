
%% 1. Process behavioral data
options.stride_sorting = 1;
options.compute_gait_parameters = 1;

options.group_parameters_by_distance = 0;
options.group_parameters_by_trial = 0;
options.group_speed_binned_parameters = 0;
options.group_speed_binned_parameters_by_trial = 0;

options.plot_adaptation = 0;
options.plot_tiedbelt = 0;

options.sort_strides_by_speed = 0;  % considers strides with absolute body speed around zero (i.e. small movement relative to the belt)
options.compute_imu_parameters = 0;
options.sort_front_hind_independently = 1; % It sets whether strides from the front limbs should be set independently from those of the hind limbs

params_list.gait = {'stance_length', 'stance_length_rel_body', 'stance_length_rel_body_and_belt', 'stance_speed', 'stance_duration', 'swing_onset_x_rel'...
                        'swing_length', 'swing_length_rel_body', 'swing_length_rel_body_and_belt', 'swing_speed', 'swing_duration', 'stance_onset_x_rel', ...
                        'cadence', 'duty_factor', 'double_support', 'stride_duration', 'step_length', 'coo_body', 'body_speed', 'body_speed_rel', 'stride_length', ...
                        'stance_phase', 'coo_stance_absolute', 'coo_swing_absolute'};

params_list.gait = {'stance_phase'};
                    
params_list.trajectory = {'swing_rel_x_traj'};
params_list.imu = {'imu_means', 'imu_stds'};

params_list.plot = {'stance_length', 'stance_duration', 'stance_speed', ...
     'swing_duration', 'swing_length', ...
     'double_support', 'step_length', 'coo_body'};

 
tracking_dir = [];
tracking_dir = 'Z:\LocomotionExperiments\RT Self-paced\LocomotorLearning\TM TRACKING FILES\20181018 - water deprived first test';

platform = get_default_widefield_rotary_treadmill_parameters(2);
[dirs] = LocomotionAnalysis(tracking_dir, platform, params_list, options);



%% 2. process imaging data

% set options
% the output file should be the parent folder of the animals
%input_path = 'Z:\LocomotionExperiments\RT Self-paced\Imaging\TM RAW FILES\voluntary locomotion\H1\S4';
%input_path = 'Z:\hugo.marques\LocomotionExperiments\RT Self-paced\Imaging\TM IMAGING FILES\voluntary locomotion\MC318\S2\concatenated.tif';
input_path = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S2';
output_path = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S2';
%input_path = output_path;

dirs = get_directory_tree_from_path(input_path);
%output_path = [output_path, filesep, dirs{end}];

options = [];

options.transform_imaging_videos = 1;
options.copy_syncronization_files = 0;
options.aggregate_imaging_files = 0;
options.remove_dark_frames = 0;
options.register_imaging_videos = 1;

options.overwrite_transform_directory = 0;
options.overwrite_registration_template = 0;
options.overwrite_registered_files = 0;
options.registered_file_suffix = '_Reg';
options.registered_final_suffix = '_Ready';
options.registration_filtered_files_suffix = '_Filtered';

options.registration_template_filename = 'registration_template.tif';
options.register_videos_in_directory_independently = 0;
options.delete_registration_filtered_files = 1;
options.delete_registered_files = 1;
options.delete_transformation_files = 1;

options = get_imaging_pipeline_options(options);

% pre-process pipeline
process_imaging_pipeline(input_path, output_path, options);

%% concatenate videos to make the analysis statistically relevant

input_concatenate = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S2';
output_concatenate = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S2';

N = 5 ; %concatenation group size

concatenate_mukamel(N,input_concatenate,output_concatenate)

%% Median filtering (Diogo) : removing the shutter noise by sacrifying temporal precision at first (3-frames sliding window)

input_median = 'MC318_26_199_M_tied_134,000_141,000_2__concatenated_trials_15_16_1_27_29.tif'; 

[~] = med_filter_calcium(input_median); 

%%
% extract rois using mukamel
%TODO: consider concatenating a few videos beforehand (see above)
options.nPCs = 400;
options.used_PCs = 1:350; % if this value is empty the algorithm prompts the user for manual selection
options.mu = 0.15;

cells_sort_file_muk = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5\MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51.tif';
cells_sort_file_reg = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5\MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51.tif'; 
cells_sort_out_folder = 'Z:\leonard.dupont\TM IMAGING PROCESSING FILES\voluntary locomotion\MC318\S4\mukamel_concatenated_mu0.15';
create_folder_path(cells_sort_out_folder);


muk = get_mukamel_rois(cells_sort_file_muk, [], options.nPCs, options.used_PCs, options.mu, []);

% convert mukamel rois to carey neuron
cn_rois = convert_mukamel_to_carey_neuron(muk.ica_segments, muk.seg_centroid);


% add bkg rois
[bkg,~] = draw_manual_rois(cells_sort_file_reg,[],[],[],[1,500]);
cn_bkg_rois = convert_manual_rois_to_carey_neuron(bkg);


% compute avg intensities for rois and background, and save data : done on
% the non-filtered tiff file so that we do not loose temp. resolution
[~, frame_avg, cn] = get_carey_neurons_mean_intensity(cells_sort_file_reg, cn_rois, [cells_sort_out_folder,'\','carey_neurons.mat']);
[~, ~, cn_bkg] = get_carey_neurons_mean_intensity(cells_sort_file_reg, cn_bkg_rois, [cells_sort_out_folder,'\','carey_neurons_bkg.mat']);

%plot the calcium traces
plot_all_traces(cn.intensity)

data = load('Z:\leonard.dupont\TM IMAGING PROCESSING FILES\voluntary locomotion\MC318\S4\mukamel_concatenated2\carey_neurons.mat');
cn = data.cn;

data_bkg = load('Z:\leonard.dupont\TM IMAGING PROCESSING FILES\voluntary locomotion\MC318\S4\mukamel_concatenated2\carey_neurons_bkg.mat');
cn_bkg = data_bkg.cn;


% Convert calcium traces to spikes - not useful if using the following
% MLSpike algorithm 
ops.fs = 30;
ops.recomputeKernel = 0;
ops.sensorTau = 0.7;
ops.estimateNeuropil = 0;
%ops.deconvType = 'L0'; 
ops.deconvType = 'OASIS';
threshold = 1000;
[~, ~, cn] = get_spikes_from_calcium_traces(cn.intensity, ops, threshold, cn, []);


%% Automatic removal of bad ROIs
% We trained an SVM to be able to automatically remove bad regions of
% interest based on shape criteria (see description of the extraction
% function). 

svmname = 'SVM_Pkj.mat'; 
embedded = load(svmname); 
TheSVM = embedded.Pkj_sorter.ClassificationSVM; 
clear embedded
criteria = extract_mask_criteria(cn); 
labels = predict(TheSVM,criteria); 

cn = remove_bad_purkinje(cn,labels); 
clear labels
clear criteria

%% Check for 'border' rois (code will come one day)
% Execute blocks sequentially

clear good_purkinje
manual_roi_sorting(cn)

good_purkinje = getglobal_purkinje;

cn = remove_bad_purkinje(cn,good_purkinje); 

%% Elliptical FISSA demixing
% March 2019 


% . . . . . . . . This if you want to include background in demixing :
% eroding might be an option if it's a large area  . . . . . . . . 
msk = cn_bkg.mask{1,1};
msk = imerode(msk,strel('disk',20));
cn_bkg.mask{1,1} = msk; 

% . . . . . . . Build ellipsoids around cell bodies (Pkj) . . . . . . . .  
nseg = 4; 
bkg = true;
bkgmask = cn_bkg.mask{1,1};
l = 7;
neuropile = ...
    building_ellipses(cn,'graphics',0,'segment_ellipse',1,'nseg',nseg,'bkg',bkg,'bkgmask',bkgmask,'l',l);


% . . . . . . . Make fluorescence traces from ellipsoid masks . . . . . . .  
inputfile = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51.tif'; 
neuropile = subR_fluorescence(neuropile,inputfile); 

plot_segcell_traces(neuropile,12) %for visualisation 

if neuropile.bkg
    nseg = nseg+1;
end


N = neuropile.n_cells;
graphics = 0;
% . . . . . . . Normalise before demixing, otherwise artefacts . . . . . . 
for cell = 1:N
    for seg = 1:(nseg+1)
        a = zero_and_max(neuropile.intensity{seg,cell}); 
        neuropile.intensity{seg,cell} = a ; 
    end
end



% . . . . . . . . NNMF occurs here with NNDSVD initialisation . . . . . . .
for cell = 1:N
    
   if rem(cell,10) == 0 || cell == N
       disp(['Cell ',num2str(cell),' out of ',num2str(N),'.'])
   end

    if graphics
        for seg = 1:nseg+1
            start = 1 + (seg-1)*3 ;
            stop = seg * 3;
            subplot(nseg+1,3,start:stop)
            plot(neuropile.intensity{seg,cell},'color',cmap(seg,:))
            box off
            start = stop + 1;
            stop = stop + 1;
        end    
    end
    
    tic
    disp('1 -- Preparing fluorescence matrix with neuropile subregions --')
    F = initFISSA_mixedF(neuropile,cell);
    toc
    disp('2 -- NNDSVD-based initialisation --')
    [W0,H0] = NNDSVD(F,nseg+1,0); 
    %initialisation of weight and separated matrices using SVD 
    

    % . . . . . . . . minimising |F - W*H| according to Frobenius . . . . .
    [W,H] = nnmf(F,nseg+1,'w0',W0,'h0',H0); 
   
    if neuropile.bkg
        maxi = max(W(nseg,:)); %then last but one segment (last = bkg)
        wc = find(W(nseg,:) == maxi);
    else
        maxi = max(W(nseg+1,:)); %then last segment
        wc = find(W(nseg+1,:) == maxi);
    end
    dmxd = H(wc,:);
    
    cn.intensity_dm(:,cell) = dmxd; 
    
    if graphics 
        figure, hold on
        axh1 = subplot(2,3,1:3);
        mixed = zero_and_max(cn.intensity(:,cell).');
        plot(mixed,'k')
        axis tight, box off
        title('mixed')
        axh2 = subplot(2,3,4:6); 
        demixed = zero_and_max(cn.intensity_dm(:,cell).');
        plot(demixed,'k')
        axis tight, box off 
        title('demixed')
        linkaxes([axh1,axh2],'xy')
        hold off
    end
 
end

%visualise, compare
compare_mx_dmx(cn,48)
plot_segcell_traces(neuropile,48)
plot_all_traces(cn.intensity)
plot_all_traces(cn.intensity_dm)

%% Clustering using k-means ++ and hierarchical merging
% February 2019


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

%% Giving a go to the MLSpike algorithm (Deneux et al. 2016) https://www.nature.com/articles/ncomms12190
% a : selecting single spike events
% we must first use a nice calcium recording vector from which we'll select
% single spike events. In this way, MLSpike will autocalibrate. 
calcium_data = cn.intensity_dmdn;  
roi = 45; %check with the plot_all_traces function
good_trace = calcium_data(:,roi);
[pks,locs] = define_single_events(good_trace);

% execute separately - GetGlobal command
sspike_events = select_single_spikes(locs);

% b : now we feed the autocalibration and deconvolute the signals 
clear deconvolution
clear par
[x,y] = size(calcium_data);


% autocalibration
calcium_norm = good_trace/mean(good_trace); %doesn't necessarily need to be normalised

pax = spk_autocalibration('par');
pax.dt = 1/30; %framerate
pax.realspikes = sspike_events*1/30; %the transients we deem as representative of sspike events (see above)
pax.mlspikepar.dographsummary = false;

[tauest, aest, sigmaest] = spk_autocalibration(calcium_norm,pax);

% deconvolution using Viterbi algorithm
par = tps_mlspikes('par');
par.a = aest;
par.dt = 1/30;
par.tau = tauest;
par.drift.parameter = 0.001;
p2 = 0.5;
p3 = 0.01;
par.pnonlin = [p2 p3] ;
par.dographsummary = false;
par.finetune.sigma = sigmaest; %to be left empty from my experience, even
%though we estimated it 


for roi=1:185
    disp(['Processing ROI ' , num2str(roi)])
    [spk, fit, drift, parest] = spk_est(calcium_data(:,roi),par); 
    cn.spiketimes{1,roi} = spk;
end

%%
% visualisation

roi=29; %choose the region (can eventually be looped if we want to check them all)
spike_times = deconvolution.(['roi_',num2str(roi)]).spiketimes; 
spk_display(par.dt,spike_times,calcium_data(:,roi))

%% c : build ISI histogram

[~] = build_ISI_histo(deconvolution,'individual',0,'bw',0.01,'rmhigh',4);

%% d : build rasterplot 

rasterplot(deconvolution)
%%

tiffimage = '/Users/leonarddupont/Desktop/M2_internship/registration_template.tif';

[spatial] = spatial_rasterplot(tiffimage,deconvolution,5); 


%% 3. Process session data

session_raw_dir = 'Z:\leonard.dupont\TM RAW FILES\voluntary locomotion\MC318\S5';
tracking_dir = 'Z:\leonard.dupont\TM TRACKED FILES\MC318\S5';
imaging_dir = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5';
session_output_dir = 'Z:\leonard.dupont\TM SESSION FILES\voluntary locomotion\';

[exp_files] = get_experimental_files_ordered_by_animal_session_and_trial(session_raw_dir, '*.tdms');
platform = get_default_widefield_rotary_treadmill_parameters(2);

concatenate_session_TDMS_data(exp_files, session_output_dir, platform);
concatenate_session_RTM_treadmill_data(session_output_dir, platform);
concatenate_session_RTM_tracking_data(exp_files, tracking_dir, session_output_dir, platform);
concatenate_session_RTM_imaging_data(exp_files, imaging_dir, session_output_dir, platform);

load('/Volumes/carey/leonard.dupont/TM SESSION FILES/voluntary locomotion/MC318/S5/imaging_data.mat')
load('/Volumes/carey/leonard.dupont/TM SESSION FILES/voluntary locomotion/MC318/S5/port_data.mat')
load('/Volumes/carey/leonard.dupont/TM SESSION FILES/voluntary locomotion/MC318/S5/tracking_data.mat')
load('/Volumes/carey/leonard.dupont/TM SESSION FILES/voluntary locomotion/MC318/S5/treadmill_data.mat')

%%
figure, hold on
plot(tm.time,tm.speedM), axis tight
calcium_data = cn.intensity_dm.';
m_ca = mean(calcium_data,1);
plot(im.time,m_ca)




figure, hold on
plot(tm.time(1:length(tm_speedM)),tm_speedMs), axis tight
plot(linspace(0,tm.time(length(tm_speedM)),length(m_ca)),m_ca)


% . . . . . . . . SPEED CATEGORIES ANALYSIS . . . . . . . . . . . . . . . . 
ncat = 4; 

tm_speedMs(tm_speedMs < 0) = 0;
maxs = max(tm_speedMs);
mins = min(tm_speedMs);

diffcat = (maxs - mins)/ncat ;
speedcats = zeros(ncat+1,1);
for i = 1:ncat+1
    speedcats(i) = mins + (i-1) * diffcat;
end

speed_labels = zeros(length(tm_speedM),1);
for j = 1:length(tm_speedM)
    cat = 1;
    while tm_speedMs(j) > speedcats(cat)
        cat = cat + 1;
    end
    speed_labels(j) = cat;
end

plot(speed_labels)

ds = 100;

cmap = jet(ncat+1);
coord = linspace(1,length(tm_speedM),length(tm_speedM));
speed_labelds = downsample(speed_labels,ds);
tm_speedMds = downsample(tm_speedMs,ds);

point = 1;
xc = [];
label = 0;
while point < length(speed_labels)
    nlabel = speed_labels(point);
    if nlabel ~= label
        xc(end+1) = point;
        label = nlabel;
    end
    point = point + 1;
end
xc(end+1) = point;


figure, hold on

% plot_dummies for legend
l1 = plot([NaN,NaN], 'color', cmap(1,:));
l2 = plot([NaN,NaN], 'color', cmap(2,:));
l3 = plot([NaN,NaN], 'color', cmap(3,:));
l4 = plot([NaN,NaN], 'color', cmap(4,:));
l5 = plot([NaN,NaN], 'color', cmap(5,:));
legend([l1, l2, l3, l4, l5], {'vlow','low','medium','fast','vfast'},'AutoUpdate','off');

for k = 1:length(xc)-1
    xplus = xc(k+1);
    xminus = xc(k);
    label = speed_labels(xplus-1);
    color = cmap(label,:);
    
    xl = [xminus, xminus, xplus, xplus];
    yl = [0 100 100 0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.4;
    a.EdgeColor = 'none';
end
        
axis tight, hold on
plot(tm_speedMs*(100/maxs),':','color','k','LineWidth',1)

imt = linspace(0,length(tm_speedMs),length(m_ca));
plot(imt,m_ca*(100/max(m_ca)),'color','k','LineWidth',0.001)
box off 
yticklabels = [];

tm_speedMds = downsample(tm_speedM,round(length(tm_speedM) / length(m_ca)));

%%
N = cn.n_cells;
synchros = zeros(N,N);

c=0;
for i = 1:N
    for j = 1:N
        c = c+1;
        percent = c/(N*N) * 100;
        if rem(percent,5) == 0 || percent == 0 
            disp([num2str(percent),'% completed.'])
        end
        
        x = cn.intensity(:,i).';
        y = cn.intensity(:,j).';
        
        x = zero_and_max(x);
        y = zero_and_max(y);
        
        [xpk,tx] = findpeaks(x,'MinPeakHeight',0.1,'MinPeakProminence',0.05);
        [ypk,ty] = findpeaks(y,'MinPeakHeight',0.1,'MinPeakProminence',0.05);

        [es,~,~] = Event_Sync(tx,ty); 
        synchros(i,j) = es;
    end
end

