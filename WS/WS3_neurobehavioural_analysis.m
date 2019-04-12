%% CareyLab - 2019 - Working Script 3
% This script uses the data extracted in WS1-2 and the speed+tracking data
% to dive into the actual analysis we are interested in.
% Still in the making. 
%
%       (A)  Tracking
%       (B)  TDMS processing 
%       (C)  Speed-based analysis
%       (D)  Point-mutual information (PMI) 
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


%% (A) - Tracking 
% Uses the Locomouse tracker to get precise behavioural information
options.stride_sorting = 1;
options.compute_gait_parameters = 1;
options.group_parameters_by_distance = 0;
options.group_parameters_by_trial = 0;
options.group_speed_binned_parameters = 0;
options.group_speed_binned_parameters_by_trial = 0;
options.plot_adaptation = 0;
options.plot_tiedbelt = 0;
options.sort_strides_by_speed = 0;  
% considers strides with absolute body speed around zero 
%(i.e. small movement relative to the belt)
options.compute_imu_parameters = 0;
options.sort_front_hind_independently = 1;
% It sets whether strides from the front limbs should be set 
% independently from those of the hind limbs

params_list.gait = ...
    {'stance_length', 'stance_length_rel_body',...
    'stance_length_rel_body_and_belt', 'stance_speed', 'stance_duration',...
    'swing_onset_x_rel','swing_length', 'swing_length_rel_body',...
    'swing_length_rel_body_and_belt', 'swing_speed', 'swing_duration',...
    'stance_onset_x_rel', 'cadence', 'duty_factor', 'double_support',...
    'stride_duration', 'step_length', 'coo_body', 'body_speed',...
    'body_speed_rel', 'stride_length', 'stance_phase',...
    'coo_stance_absolute', 'coo_swing_absolute'};

params_list.gait = {'stance_phase'};
                    
params_list.trajectory = {'swing_rel_x_traj'};
params_list.imu = {'imu_means', 'imu_stds'};

params_list.plot = {'stance_length', 'stance_duration', 'stance_speed', ...
     'swing_duration', 'swing_length', ...
     'double_support', 'step_length', 'coo_body'};

 
tracking_dir = [];
tracking_dir = 'Z:\LocomotionExperiments\RT Self-paced\LocomotorLearning\TM TRACKING FILES\20181018 - water deprived first test';

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
platform = get_default_widefield_rotary_treadmill_parameters(2);
[dirs] = LocomotionAnalysis(tracking_dir, platform, params_list, options);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% (B) Processing session TDMS files
% TDMS files are metadata files that are created all along the session. TTL
% inputs are written in distinct files to have sampling-frequency
% compatible timestamps across all devices (microscope, behaviour camera,
% belt speed, etc.). 
% Outputs are im, tm, tracking and port_data structs. The 'time' fields can
% be used to have a common timeframe. 


session_raw_dir = 'Z:\leonard.dupont\TM RAW FILES\voluntary locomotion\MC318\S5';
tracking_dir = 'Z:\leonard.dupont\TM TRACKED FILES\MC318\S5';
imaging_dir = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5';
session_output_dir = 'Z:\leonard.dupont\TM SESSION FILES\voluntary locomotion\';

[exp_files] = get_experimental_files_ordered_by_animal_session_and_trial(session_raw_dir, '*.tdms');
platform = get_default_widefield_rotary_treadmill_parameters(2);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
concatenate_session_TDMS_data(exp_files, session_output_dir, platform);
concatenate_session_RTM_treadmill_data(session_output_dir, platform);
concatenate_session_RTM_tracking_data(exp_files, tracking_dir, session_output_dir, platform);
concatenate_session_RTM_imaging_data(exp_files, imaging_dir, session_output_dir, platform);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% (C) Analysing belt speed
% Hereunder is a rough draft to co-analyse speed variations on the RTM 
% and calcium transients. 
% (C1) First, we take a look at the plot over the whole
% session. Since imaging trials might have been concatenated and separately
% analyse afterwards, we need to concatenate speed trials accordingly. We
% also smooth the speed for readability and analysis purposes. 
% (C2) Then, we define speed categories that we will be able to use later
% one when correlating fluorescence to behaviour. 


% C1 - speed data preprocessing 
calcium_data = cn.intensity.';
m_ca = mean(calcium_data,1);
m_ca = zero_and_max(m_ca);

% Initial plot : over the whole session 
figure, hold on
plot(tm.time,tm.speedM), axis tight
plot(im.time,m_ca)

Nkk = 4; % SAME AS CONCATENATE MUKAMEL 
inputpath = '';
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
speedkk = imbased_konkat(im,tm,inputpath,Nkk); 
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

if isstruct(speedkk)
    speedkk_struct = speedkk;
    speedkk = speedkk_struct.kk(1); %take the first one 
end

tm_speedMs = smoothdata(speedkk,'gaussian',10000); %smoothed version
figure, plot(tm_speedMs)

spdtime = linspace(1,tm.time(length(tm_speedMs)),length(tm_speedMs)); 
imtime = linspace(1,tm.time(length(tm_speedMs)),length(m_ca));

rescale = max(tm_speedMs);

figure, hold on
plot(spdtime,tm_speedMs), axis tight
plot(imtime,m_ca*rescale)


% C2 - analysing speed categories 
ncat = 2; 
% - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
speedcats = define_speedcats(tm_speedMs,ncat);
% - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
speedlabels = assign_speedlabels(tm_speedMs,speedcats);
ncat = max(speedlabels); % speed == 0 always adds a category, check if so
% - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
xc = get_speedregime_boundaries(speedlabels);
% - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

plot_all_traces_w_speed(cn.intensity,xc,speedlabels,tm_speedMs)

%% (D) - Relating fluorescence to speed
% What we now want is to see how changes in fluorescence (cor)relate to
% changes in belt speed, id est in the animal's behavioural state (coarse
% window into it, not tracking yet). 
% For a start we do a mutual-information based analysis. 

Tim = length(m_ca);
Ttm = length(tm_speedMs);
mergeF = Ttm/Tim; 

binboundaries = zeros(k,1); 
for k = 1:Tim
    binboundaries(k) = floor((k-1)*mergeF + 1);
end

aside = Ttm - binboundaries(end);
binboundaries(end+1) = binboundaries(end)+aside;

figure, hold on
Nclust = length(activity_cluster.clusterregions);
for j = 1:Nclust
    regions = activity_cluster.clusterregions{:,j};
    N = length(regions);
    allclusttraces = zeros(length(cn.intensity(:,1)),N);
    for roi = 1:N
        allclusttraces(:,roi) = zero_and_max(cn.intensity(:,regions(roi)).').';
    end
    mclusttrace = mean(allclusttraces,2);
    mclusttrace = mclusttrace.';
    %for roi = 1:N
        %trace = cn.intensity(:,roi).';
        %trace = zero_and_max(trace);
        fluocats = define_fluocats(mclusttrace,10);
        fluolabels = assign_fluolabels(mclusttrace,fluocats);

        Nimbins = max(fluolabels);
        Ntmbins = max(speedlabels);
        MImat = zeros(Nimbins,Ntmbins);

        for k = 1:Tim
            start = binboundaries(k);
            stop = binboundaries(k+1);
            tmlabel = round(mean(speedlabels(start:stop)));
            imlabel = fluolabels(k);
            MImat(imlabel,tmlabel) = MImat(imlabel,tmlabel) + 1;
        end

        MImat = MImat/Tim; 

        for imlbl = 1:Nimbins
            Pim = length(find(fluolabels == imlbl))/ Tim;
            for tmlbl = 1:Ntmbins
                Ptm = length(find(speedlabels == tmlbl))/ Ttm;
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                MImat(imlbl,tmlbl) =  log(MImat(imlbl,tmlbl)/ (Pim*Ptm));
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            end
        end





    subplot(2,3,j),imagesc(MImat)
    title(['Cluster ',num2str(j)])
    ylabel('Fluorescence bins')
    xlabel('Speedbins')
    if rem(j,Nclust/2) == 0
        colorbar
    end
end
suptitle('Point mutual information in all clusters')


%% Consider saving the workspace?

wkspacename = ['WS1_workspace',WSanalysis.dd,WSanalysis.mouse,WSanalysis.session]; 
save(wkspacename)