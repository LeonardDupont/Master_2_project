%% CareyLab - 2019 - Working Script 1
% This script uses 1P fluorescence data recorded in vivo and processes it
% with the hereunder-mentioned steps:
%
%       (A)  Image registration to queen frame (registration template)
%       (B)  Trial concatenation in groups of (Nkk) videos in order to use
%            statistically relevant samples for segmentation.
%       (C)  Segmentation (Mukamel et al. 2009), coupled PCA-stICA.
%       (D)  Background extraction, structure (carey neuron) creation.
%       (E)  ROI sorting using either a trained SVM or manual sorting, or
%            both. 
%      (PL)  Plotting tools 
%
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


%% (A) Image registration to queen frame


input_path = ...
    'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S2';
output_path = ...
    'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S2';

dirs = get_directory_tree_from_path(input_path);

% - - - - - - - - - - - registration options - - - - - - - - - - - - - - - 
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
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
process_imaging_pipeline(input_path, output_path, options);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% (B) Trial concatenation (kk)


input_kk = ...
    'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S2';
output_kk = ...
    'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S2';
Nkk = 5;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
concatenate_mukamel(Nkk,input_kk,output_kk)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% (C) Segmentation, Mukamel

cells_sort_file_muk = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5\MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51.tif';
cells_sort_file_reg = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5\MC318_26_218_F_tied_138,000_137,000_5__concatenated_trials_1_25_27_51.tif'; 
cells_sort_out_folder = 'Z:\leonard.dupont\TM IMAGING PROCESSING FILES\voluntary locomotion\MC318\S4\mukamel_concatenated_mu0.15';
create_folder_path(cells_sort_out_folder);

nPCs = 400;
used_PCs = 1:350; 
mu = 0.15; %we take a mainly spatial ICA, as recommended 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
muk = get_mukamel_rois(cells_sort_file_muk, [], nPCs, used_PCs, mu, []);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% (D) Creating carey neuron structures, creating background ROI

% D1 - regions of interest (neurons)
cn_rois = convert_mukamel_to_carey_neuron(muk.ica_segments, muk.seg_centroid);
[~, frame_avg, cn] = get_carey_neurons_mean_intensity(cells_sort_file_reg, cn_rois, [cells_sort_out_folder,'\','carey_neurons.mat']);
data = load('Z:\leonard.dupont\TM IMAGING PROCESSING FILES\voluntary locomotion\MC318\S4\mukamel_concatenated2\carey_neurons.mat');
cn = data.cn;

% D2 - background if needed
    %creating
[bkg,~] = draw_manual_rois(cells_sort_file_reg,[],[],[],[1,500]);
cn_bkg_rois = convert_manual_rois_to_carey_neuron(bkg);

    %saving
[~, ~, cn_bkg] = get_carey_neurons_mean_intensity(cells_sort_file_reg, cn_bkg_rois, [cells_sort_out_folder,'\','carey_neurons_bkg.mat']);
data_bkg = load('Z:\leonard.dupont\TM IMAGING PROCESSING FILES\voluntary locomotion\MC318\S4\mukamel_concatenated2\carey_neurons_bkg.mat');
cn_bkg = data_bkg.cn;


%% (E) Selecting rois to keep the relevant ones

% E1 - Automatic with support vector machine 
% Saving results in an independent struct is recommended at first
svmname = 'SVM_Pkj.mat'; 
embedded = load(svmname); 
TheSVM = embedded.Pkj_sorter.ClassificationSVM; 
clear embedded
criteria = extract_mask_criteria(cn); 
labels = predict(TheSVM,criteria); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - 
cn = remove_bad_purkinje(cn,labels); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - 
clear labels
clear criteria

% E2 - Manual removal of rois
% Use is advised even after automatic removal to check for border ROIs
% Execute all steps independently
clear good_purkinje
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
manual_roi_sorting(cn)
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
good_purkinje = getglobal_purkinje;
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
cn = remove_bad_purkinje(cn,good_purkinje); 
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -


%% (PL) Plotting tools 

%visualise masks (before and after roi sorting e.g.)
purkinje_artscape(cn.mask)
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
artfigname = ['Pkj_artscape',WSanalysis.dd,WSanalysis.mouse,WSanalysis.session];
print(artfigname,'-dpdf','-r500')



%plot all calcium traces in a single window
plot_all_traces(cn.intensity)
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
alltracesname = ['All_traces',WSanalysis.dd,WSanalysis.mouse,WSanalysis.session]; 
print(alltracesname,'-dpdf','-r500') 


%% Consider saving the workspace?

wkspacename = ['WS1_workspace',WSanalysis.dd,WSanalysis.mouse,WSanalysis.session]; 
save(wkspacename)








