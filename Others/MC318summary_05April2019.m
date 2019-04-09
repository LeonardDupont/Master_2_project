%% DATA ANALYSIS SUMMARY - April 2019
% This whole script was written as a summary script to make general figures
% being part of a progress report. All datafiles come from MC318, a male
% mice from the L7-Cre line that was injected with
% AAV1/2.CAG.GCaMP6f(floxed) in the right paravermal cortex. 

% leonard.dupont@ens.fr
%
%% Data loading
data = load('/Users/leonarddupont/Desktop/M2_internship/Data/S2/carey_neurons.mat');
cns2 = data.cn;

data = load('/Users/leonarddupont/Desktop/M2_internship/Data/S4/carey_neurons.mat');
cns4 = data.cn;

data = load('/Users/leonarddupont/Desktop/M2_internship/Data/S5/carey_neurons.mat');
cns5 = data.cn;

clear data 
nS = 3; %number of sessions 

%% Plotting the masks for all sessions, before and after selection
% E1 - Automatic with support vector machine 
% Saving results in an independent struct is recommended at first
svmname = 'SVM_Pkj.mat'; 
embedded = load(svmname); 
TheSVM = embedded.Pkj_sorter.ClassificationSVM; 
clear embedded
clear svmname
criteria = extract_mask_criteria(cns2s); 
labels = predict(TheSVM,criteria); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - 
cns2s = remove_bad_purkinje(cns2s,labels); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - 
clear labels
clear criteria

% E2 - Manual removal of rois
% Use is advised even after automatic removal to check for border ROIs
% Execute all steps independently
clear good_purkinje
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
manual_roi_sorting(cns5s)
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
good_purkinje = getglobal_purkinje;
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
cns4s = remove_bad_purkinje(cns4,good_purkinje); 
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

%% PLOTS THEMSELVES 
figure, hold on
subplot(2,3,1), purkinje_artscape(cns2.mask)
title('Session 2')
subplot(2,3,2), purkinje_artscape(cns4.mask)
title('Session 4')
subplot(2,3,3), purkinje_artscape(cns5.mask)
title('Session 5') 
subplot(2,3,4), purkinje_artscape(cns2s.mask)
subplot(2,3,5), purkinje_artscape(cns4s.mask)
subplot(2,3,6), purkinje_artscape(cns5s.mask)

%%
nrois = 3; 
S2 = [21,23,5];
S4 = [13,65,16];
S5 = [5,78,16];

figure, hold on
subplot(1,3,1),purkinje_artscape(cns2s.mask,S2), title('Session 2')
subplot(1,3,2),purkinje_artscape(cns4s.mask,S4), title('Session 4')
subplot(1,3,3),purkinje_artscape(cns5s.mask,S5), title('Session 5') 


%% Manual spike count
% - - - - - - - - - - - - - S2 manual spikes  - - - - - - - - - - - - - - -
clear manual_spikes 
manually_count_spikes(cns2s.intensity(:,S2(1)))
spiketimes = getGlobal_manualspikes; 
timings = zeros(1,length(cns2s.intensity(:,S2(1))));
timings(spiketimes) = 1;
cns2s.mspikes(1,:) = timings; 

clear manual_spikes 
manually_count_spikes(cns2s.intensity(:,S2(2)))
spiketimes = getGlobal_manualspikes; 
timings = zeros(1,length(cns2s.intensity(:,S2(1))));
timings(spiketimes) = 1;
cns2s.mspikes(2,:) = timings; 

clear manual_spikes 
manually_count_spikes(cns2s.intensity(:,S2(3)))
spiketimes = getGlobal_manualspikes; 
timings = zeros(1,length(cns2s.intensity(:,S2(1))));
timings(spiketimes) = 1;
cns2s.mspikes(3,:) = timings; 

% - - - - - - - - - S4 manual spikes  - - - - - - - - - - - - - - - - - - -
clear manual_spikes 
manually_count_spikes(cns4s.intensity(:,S4(1)))
spiketimes = getGlobal_manualspikes; 
timings = zeros(1,length(cns4s.intensity(:,S4(1))));
timings(spiketimes) = 1;
cns4s.mspikes(1,:) = timings; 

clear manual_spikes 
manually_count_spikes(cns4s.intensity(:,S4(2)))
spiketimes = getGlobal_manualspikes; 
timings = zeros(1,length(cns4s.intensity(:,S4(1))));
timings(spiketimes) = 1;
cns4s.mspikes(2,:) = timings; 

clear manual_spikes 
manually_count_spikes(cns4s.intensity(:,S4(3)))
spiketimes = getGlobal_manualspikes; 
timings = zeros(1,length(cns4s.intensity(:,S4(1))));
timings(spiketimes) = 1;
cns4s.mspikes(3,:) = timings; 

% - - - - - - - - - S5 manual spikes  - - - - - - - - - - - - - - - - - - -
clear manual_spikes 
manually_count_spikes(cns5s.intensity(:,S5(1)))
spiketimes = getGlobal_manualspikes; 
timings = zeros(1,length(cns5s.intensity(:,S5(1))));
timings(spiketimes) = 1;
cns5s.mspikes(1,:) = timings; 

clear manual_spikes 
manually_count_spikes(cns5s.intensity(:,S5(2)))
spiketimes = getGlobal_manualspikes; 
timings = zeros(1,length(cns5s.intensity(:,S5(1))));
timings(spiketimes) = 1;
cns5s.mspikes(2,:) = timings; 

clear manual_spikes 
manually_count_spikes(cns5s.intensity(:,S5(3)))
spiketimes = getGlobal_manualspikes; 
timings = zeros(1,length(cns5s.intensity(:,S5(1))));
timings(spiketimes) = 1;
cns5s.mspikes(3,:) = timings; 

%% Automatic spike count

ops.fs = 30;
ops.recomputeKernel =0;
ops.sensorTau = 0.7;
ops.estimateNeuropil = 0;
ops.deconvType = 'L0'; 
%ops.deconvType = 'OASIS';
threshold = se_amp;
[~, ~, ~, cns2s] = get_spikes_from_calcium_traces(cns2s.intensity, ops, threshold, cns2s, []);
[~, ~, ~, cns4s] = get_spikes_from_calcium_traces(cns4s.intensity, ops, threshold, cns4s, []);
[~, ~, ~, cns5s] = get_spikes_from_calcium_traces(cns5s.intensity, ops, threshold, cns5s, []);

%%
cmap = parula(5);
sizes2 = size(cns2s.mask{1,1});
sizes4 = size(cns4s.mask{1,1});
sizes5 = size(cns5s.mask{1,1});
pt_max = 3497;
time = linspace(1,pt_max/30,pt_max); 
%%

for i = 1:nrois
    figure('Renderer', 'painters', 'Position', [500 500 900 600])
    subplot(4,3,1)
    
        S = sizes2;
        white = cat(3,zeros(S),zeros(S),zeros(S)); 
        h = imshow(white); hold on 
        art = S2(i);  
        I = cns2s.mask{1,art};
        c = cmap(i,:);
        full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
        h = imshow(full); hold on 
        set(h, 'AlphaData', I) , hold off 
        title('Session 2')
    
    subplot(4,3,2)
    
        S = sizes4;
        white = cat(3,zeros(S),zeros(S),zeros(S)); 
        h = imshow(white); hold on 
        art = S4(i);  
        I = cns4s.mask{1,art};
        c = cmap(i,:);
        full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
        h = imshow(full); hold on 
        set(h, 'AlphaData', I) , hold off 
        title('Session 4')
    
    
    subplot(4,3,3)
    
        S = sizes5;
        white = cat(3,zeros(S),zeros(S),zeros(S)); 
        h = imshow(white); hold on 
        art = S5(i);  
        I = cns5s.mask{1,art};
        c = cmap(i,:);
        full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
        h = imshow(full); hold on 
        set(h, 'AlphaData', I) , hold off 
        title('Session 5') 
        
    xtstart = 60;
    xtstop = 80;
    x = [xtstart,xtstart,xtstop,xtstop];
    y = [0, 1, 1, 0];
    subplot(4,3,4)
        plot(time,zero_and_max(cns2s.intensity(1:pt_max,S2(i)).'),'color',cmap(i,:)), hold on
        u = fill(x,y,cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        xlabel('Time (s)')
        ylabel('\Delta F / F')
        axis tight
        
    subplot(4,3,5)
        plot(time,zero_and_max(cns4s.intensity(1:pt_max,S4(i)).'),'color',cmap(i,:)), hold on
        u = fill(x,y,cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        xlabel('Time (s)')
        axis tight
        
    subplot(4,3,6)
        plot(time,zero_and_max(cns5s.intensity(1:pt_max,S5(i)).'), 'color',cmap(i,:)), hold on 
        u = fill(x,y,cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        xlabel('Time (s)')
        axis tight
        
    tstrt = round(60*30);
    tstp = round(80*30); 
    time2 = linspace(60,80,tstp-tstrt+1);
    subplot(4,3,7)
        data = zero_and_max(cns2s.intensity.').'; 
        plot(time2,data(tstrt:tstp,S2(i)), 'color',cmap(i,:)), hold on
        axis tight
        u = fill([60 60 80 80],[0, max(data(tstrt:tstp,S2(i)))+0.1, max(data(tstrt:tstp,S2(i)))+0.1, 0],cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        xlabel('Time (s)')
        ylabel('\Delta F / F')
    subplot(4,3,8)
        data = zero_and_max(cns4s.intensity.').';
        plot(time2,data(tstrt:tstp,S4(i)), 'color',cmap(i,:)), hold on
        axis tight
        u = fill([60 60 80 80],[0, max(data(tstrt:tstp,S4(i)))+0.1, max(data(tstrt:tstp,S4(i)))+0.1, 0],cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        xlabel('Time (s)')
    subplot(4,3,9)
        data = zero_and_max(cns5s.intensity.').';
        plot(time2,data(tstrt:tstp,S5(i)), 'color',cmap(i,:)), hold on
        axis tight
        u = fill([60 60 80 80],[0, max(data(tstrt:tstp,S5(i)))+0.1, max(data(tstrt:tstp,S5(i)))+0.1, 0],cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        xlabel('Time (s)')
    
    subplot(4,3,10)
        ISIs = simple_ISIs(cns2s.mspikes(i,:));
        h = histogram(ISIs,'BinWidth',0.5,'Normalization','probability'); hold on
        wfreqs = slidingwindow_freq(cns2s.mspikes(i,:));
        mm = sum(wfreqs)/length(wfreqs~=0);
        ttt = ['Mean = ',num2str(mm)];
        text(max(wfreqs)/2,0.2,ttt)
        h.FaceColor = cmap(i,:);
        h.FaceAlpha = 0.8;
        ISIs2 = simple_ISIs(cns5s.spikes(:,S5(i)));
        g = histogram(ISIs2,'BinWidth',0.5,'Normalization','probability');
        g.FaceColor = cmap(i,:);
        g.FaceAlpha = 0.3; 
        xlabel('Hz bin')
        ylabel('Probability')
    subplot(4,3,11)
        ISIs = simple_ISIs(cns4s.mspikes(i,:));
        h = histogram(ISIs,'BinWidth',0.5,'Normalization','probability'); hold on
        wfreqs = slidingwindow_freq(cns4s.mspikes(i,:));
        mm = sum(wfreqs)/length(wfreqs~=0);
        ttt = ['Mean = ',num2str(mm)];
        text(max(wfreqs)/2,0.2,ttt)
        h.FaceColor = cmap(i,:);
        h.FaceAlpha = 0.8;
        ISIs2 = simple_ISIs(cns5s.spikes(:,S5(i)));
        g = histogram(ISIs2,'BinWidth',0.5,'Normalization','probability');
        g.FaceColor = cmap(i,:);
        g.FaceAlpha = 0.3; 
        xlabel('Hz bin')
    subplot(4,3,12)
        ISIs = simple_ISIs(cns5s.mspikes(i,:));
        h = histogram(ISIs,'BinWidth',0.5,'Normalization','probability'); hold on
        wfreqs = slidingwindow_freq(cns5s.mspikes(i,:));
        mm = sum(wfreqs)/length(wfreqs~=0);
        ttt = ['Mean = ',num2str(mm)];
        text(max(wfreqs)/2,0.2,ttt)
        h.FaceColor = cmap(i,:);
        h.FaceAlpha = 0.8; hold on
        ISIs2 = simple_ISIs(cns5s.spikes(:,S5(i)));
        g = histogram(ISIs2,'BinWidth',0.5,'Normalization','probability');
        g.FaceColor = cmap(i,:);
        g.FaceAlpha = 0.3; 
        xlabel('Hz bin')
end

%% Plot all traces
plot_all_traces(cns2s.intensity), hold on, title('Session 2'), hold off 
plot_all_traces(cns4s.intensity), hold on, title('Session 4'), hold off 
plot_all_traces(cns5s.intensity), hold on, title('Session 5'), hold off 

%% Clustering

clear opts
opts2.framepath = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/registration_templateS2.tif';
opts4.framepath = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/registration_templateS4.tif' ;
opts5.framepath = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/registration_template copy.tif';


clear grphcs
grphcs.dendrogram = 1;
grphcs.orgdistMAT = 1;
grphcs.orgchanceMAT = 1;
grphcs.porgchanceMAT = 1;
grphcs.clusteredlandscape = 1;

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
synchresultS2 = synchrony_clustering(cns2s,opts2,grphcs); 
synchresultS4 = synchrony_clustering(cns4s,opts4,grphcs); 
synchresultS5 = synchrony_clustering(cns5s,opts5,grphcs); 
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


%% 
figure
plot(zero_and_max(cns2s.intensity(:,S2(2)).')), hold on
plot(zero_and_max(cns2s.intensity(:,S2(3)).'))

%% 
%% (A) Non-negative matrix factorisation to demix 1P calcium signals
% All explanations are found within the used functions. 

% A1 - This if you want to include background in demixing :
% eroding might be an option if it's a large area  
msk = cns2_bkg.mask{1,1};
figure, imshow(msk)
msk = imerode(msk,strel('disk',20));
cn_bkg.mask{1,1} = msk; 


% A2 - Build ellipsoids around cell bodies - create be options 
nseg = 8; 
bkg = true;
bkgmask = cn_bkg.mask{1,1};
l = 5;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
neuropile = ...
    building_ellipses(cns2s,'graphics',0,'segment_ellipse',1,'nseg',nseg,'bkg',bkg,'bkgmask',bkgmask,'l',l);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% A3 - Extract fluorescence values corresponding to the elliptic masks
inputfile = '/Volumes/carey/hugo.marques/LocomotionExperiments/RT Self-paced/Imaging/TM IMAGING FILES/voluntary locomotion/MC318/S2/concatenated.tif';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
neuropiles2 = subR_fluorescence(neuropile,inputfile); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% A4 - Demixing using NNMF
normalise = 1;
cns2s = FISSAdemix_from_ellipses(neuropiles2,cns2s,normalise); 




%% LOAD BEHAVIOURAL DATA
% Replicate for each session
load('port_data.mat')
port_datas2 = port_data ; 
clear port_data

load('treadmill_data.mat')
tm2 = tm; 
clear tm 

load('imaging_data.mat')
im2 = im;
clear im 

%% SPEED-BASED ANALYSIS

% C1 - speed data preprocessing 
m_cas4 = mean(cns4s.intensity.',1);
m_cas4 = zero_and_max(m_cas4);

% Initial plot : over the whole session 
figure, hold on
plot(tm4.time,tm4.speedM), axis tight
plot(im4.time,m_cas4)

Nkk = 5;
inputpath = '/Volumes/carey/hugo.marques/LocomotionExperiments/RT Self-paced/Imaging/TM IMAGING FILES/voluntary locomotion/MC318/S4';
speedkk4 = imbased_konkat(im4,tm4,inputpath,Nkk);


%%

%%



tm_speedMs = speedkks4;
tm_speedMs(tm_speedMs < 0) = 0;
tm_speedMs = smoothdata(speedkks4,'gaussian',10000); %smoothed version
figure, plot(tm_speedMs)

spdtime = linspace(1,tm.time(length(tm_speedMs)),length(tm_speedMs)); 
imtime = linspace(1,tm.time(length(tm_speedMs)),length(m_cas4));

rescale = max(tm_speedMs);

figure, hold on
plot(spdtime,tm_speedMs), axis tight
plot(imtime,m_cas4*rescale)


%%
ncat = 2; 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
speedcats = define_speedcats(tm_speedMs,ncat);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
speedlabels = assign_speedlabels(tm_speedMs,speedcats);
ncat = max(speedlabels); % speed == 0 always adds a category, check if so
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
xc = get_speedregime_boundaries(speedlabels);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot_all_traces_w_speed(cns4s.intensity,xc,speedlabels,tm_speedMs)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -