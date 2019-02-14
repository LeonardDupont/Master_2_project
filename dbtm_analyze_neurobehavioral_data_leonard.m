
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
input_path = 'Z:\leonard.dupont\TM RAW FILES\voluntary locomotion\MC318\S5';
output_path = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion';
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

input_concatenate = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5';
output_concatenate = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5';

N = 5 ; %concatenation group size

concatenate_mukamel(N,input_concatenate,output_concatenate)

%% Median filtering (Diogo) : removing the shutter noise by sacrifying temporal precision at first (3-frames sliding window)

input_median = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5\MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51.tif'; 

[~] = med_filter_calcium(input_median); 

%%
% extract rois using mukamel
%TODO: consider concatenating a few videos beforehand (see above)
options.nPCs = 400;
options.used_PCs = 1:350; % if this value is empty the algorithm prompts the user for manual selection
options.mu = 0.5;

cells_sort_file_muk = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5\MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51_MEDFILT.tif';
cells_sort_file_reg = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5\MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51.tif'; 
cells_sort_out_folder = 'Z:\leonard.dupont\TM IMAGING PROCESSING FILES\voluntary locomotion\MC318\S5\mukamel_concatenated2';
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
plot_all_traces(cn)

data = load('Z:\leonard.dupont\TM IMAGING PROCESSING FILES\voluntary locomotion\MC318\S5\mukamel_concatenated2\carey_neurons.mat');
cn = data.cn;

data_bkg = load('Z:\leonard.dupont\TM IMAGING PROCESSING FILES\voluntary locomotion\MC318\S5\mukamel_concatenated2\carey_neurons_bkg.mat');
cn_bkg = data_bkg.cn;


% Convert calcium traces to spikes - not useful if using the following
% MLSpike algorithm 
ops.fs = 30;
ops.recomputeKernel = 0;
ops.sensorTau = 0.7;
ops.estimateNeuropil = 1;
ops.deconvType = 'L0'; 
%ops.deconvType = 'OASIS';
threshold = 1000;
[~, ~, cn] = get_spikes_from_calcium_traces(cn.intensity, ops, threshold, cn, []);


%% Process mukamel ROIs with CNFM to demix and then deconvolve
% DEMIXING TO BE DONE
%% Giving a go to the MLSpike algorithm (Deneux et al. 2016) https://www.nature.com/articles/ncomms12190
%% a : selecting single spike events
% we must first use a nice calcium recording vector from which we'll select
% single spike events. In this way, MLSpike will autocalibrate. 
calcium_data = cn.intensity; 
calcium_data = calcium_data(:,1:142); 
roi = 13; %check with the plot_all_traces function
good_trace = calcium_data(:,roi);
[pks,locs] = define_single_events(good_trace);

%%
% execute separately - GetGlobal command
sspike_events = select_single_spikes(locs);

%% b : now we feed the autocalibration and deconvolute the signals 
clear deconvolution
clear par
[x,y] = size(cn);


% autocalibration
calcium_norm = good_trace/mean(good_trace); %doesn't necessarily need to be normalised

pax = spk_autocalibration('par');
pax.dt = 1/30; %framerate
pax.amin = 0.003; 
pax.amax = 1;
pax.taumin = 0.1;
pax.taumax = 1.8;
pax.saturation = 7e-4; %GCaMP6f
pax.realspikes = sspike_events*pax.dt; %the transients we deem as representative of sspike events (see above)
pax.mlspikepar.dographsummary = false;

[tauest, aest, sigmaest] = spk_autocalibration(calcium_norm,pax);



% deconvolution using Viterbi algorithm
par = tps_mlspikes('par');
par.a = aest;
par.dt = 1/30;
par.tau = tauest;
par.drift.parameter = 0.0001;
par.dographsummary = false;
%par.finetune.sigma = sigmaest; to be left empty from my experience, even
%though we estimated it 


for roi=1:y 
    [spk, fit, drift, parest] = spk_est(calcium_data(:,roi),par); 
    deconvolution.(['roi_',num2str(roi)]).spiketimes = spk;  %we store everything in a big struct with rois as subfields
    deconvolution.(['roi_',num2str(roi)]).fit = fit;
    %deconvolution.(['roi_',num2str(roi)]).drift = drift; 
    deconvolution.(['roi_',num2str(roi)]).sigma = parest.finetune.sigma; 
end


% visualisation

roi = 3; %choose the region (can eventually be looped if we want to check them all)
spike_times = deconvolution.(['roi_',num2str(roi)]).spiketimes; 
spk_display(par.dt,spike_times,calcium_data(:,roi))

%% c : build ISI histogram

[~] = build_ISI_histo(deconvolution,'individual',0,'bw',0.2);


%% 3. Process session data

session_raw_dir = 'Z:\hugo.marques\LocomotionExperiments\RT Self-paced\Imaging\TM RAW FILES\voluntary locomotion\MC318\S4';
tracking_dir = 'Z:\hugo.marques\LocomotionExperiments\RT Self-paced\Imaging\TM TRACKING FILES\voluntary locomotion\MC318\S4';
imaging_dir = 'Z:\hugo.marques\LocomotionExperiments\RT Self-paced\Imaging\TM IMAGING FILES\voluntary locomotion\MC318\S4';
session_output_dir = 'Z:\hugo.marques\LocomotionExperiments\RT Self-paced\Imaging\TM SESSION FILES\voluntary locomotion\';

[exp_files] = get_experimental_files_ordered_by_animal_session_and_trial(session_raw_dir, '*.tdms');
platform = get_default_widefield_rotary_treadmill_parameters(2);

concatenate_session_TDMS_data(exp_files, session_output_dir, platform);
concatenate_session_RTM_treadmill_data(session_output_dir, platform);
concatenate_session_RTM_tracking_data(exp_files, tracking_dir, session_output_dir, platform);
concatenate_session_RTM_imaging_data(exp_files, imaging_dir, session_output_dir, platform);





