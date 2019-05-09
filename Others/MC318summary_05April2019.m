%% DATA ANALYSIS SUMMARY - April 2019
% This whole script was written as a summary script to make general figures
% being part of a progress report. All datafiles come from MC318, a male
% mice from the L7-Cre line that was injected with
% AAV1/2.CAG.GCaMP6f(floxed) in the right paravermal cortex. 

% leonard.dupont@ens.fr
%
%% Data loading
data = load('/Users/leonarddupont/Desktop/M2_internship/Data/S2/second_try/carey_neurons.mat');
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
criteria = extract_mask_criteria(cns2); 
labels = predict(TheSVM,criteria); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - 
cns2 = remove_bad_purkinje(cns2,labels); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - 
clear labels
clear criteria

% E2 - Manual removal of rois
% Use is advised even after automatic removal to check for border ROIs
% Execute all steps independently
clear good_purkinje
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
manual_roi_sorting(cns4)
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
good_purkinje = getglobal_purkinje;
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
cns4 = remove_bad_purkinje(cns4,good_purkinje); 
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
%S2 = [16,53,2];
S2 = [4,52,2];
S4 = [13,65,16];
S5 = [5,78,16];

figure, hold on
subplot(1,3,1),purkinje_artscape(cns2s.mask,S2), title('Session 2')
subplot(1,3,2),purkinje_artscape(cns4s.mask,S4), title('Session 4')
subplot(1,3,3),purkinje_artscape(cns5s.mask,S5), title('Session 5') 


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
pt_max = min([length(cns5s.intensity),length(cns4s.intensity),length(cns2s.intensity)]); 
time = linspace(1,pt_max/30,pt_max); 
%% This script bit plots a mask on top, then the fluorescence trace with the deconvolution and a zoom on it under
% It does so for all 3 selected cells in all 3 sessions. 

offset = 1;
for i = 1:nrois
    figure('Renderer', 'painters', 'Position', [500 500 900 600])
    subplot(4,3,1)
    
        S = sizes2;
        white = categ(3,zeros(S),zeros(S),zeros(S)); 
        h = imshow(white); hold on 
        art = S2(i);  
        I = cns2s.mask{1,art};
        c = cmap(i,:);
        full = categ(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
        h = imshow(full); hold on 
        set(h, 'AlphaData', I) , hold off 
        title('Session 2')
    
    subplot(4,3,2)
    
        S = sizes4;
        white = categ(3,zeros(S),zeros(S),zeros(S)); 
        h = imshow(white); hold on 
        art = S4(i);  
        I = cns4s.mask{1,art};
        c = cmap(i,:);
        full = categ(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
        h = imshow(full); hold on 
        set(h, 'AlphaData', I) , hold off 
        title('Session 4')
    
    
    subplot(4,3,3)
    
        S = sizes5;
        white = categ(3,zeros(S),zeros(S),zeros(S)); 
        h = imshow(white); hold on 
        art = S5(i);  
        I = cns5s.mask{1,art};
        c = cmap(i,:);
        full = categ(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
        h = imshow(full); hold on 
        set(h, 'AlphaData', I) , hold off 
        title('Session 5') 
        
    xtstart = 2000;
    xtstop = 2500;
    x = [xtstart,xtstart,xtstop,xtstop];
    y = [-0.2, 1, 1, -0.2];
    timeticks = linspace(1,length(cns2s.intensity(1:pt_max,1)),8);
    xlabels = cell(8,1);
    for k = 1:8
        xlabels(k) = num2cell(round(timeticks(k)/30));
    end
    subplot(4,3,4)
        plot(zero_and_max(cns2s.intensity(1:pt_max,S2(i)).'),'color',cmap(i,:)), hold on
        u = fill(x,y,cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none'; hold on
        ylim([-0.2,max(zero_and_max(cns2s.intensity(1:pt_max,S2(i)).'))+0.1])
        spikes = cns2s.spikes(1:pt_max,S2(i)).'; 
        spktimes = find(spikes==1); 
        Nspk = length(spktimes);
        for ii = 1:Nspk
            xs = [spktimes(ii),spktimes(ii)];
            ys = [-0.15,-0.05];
            plot(xs,ys,'color','k'), hold on
        end
        xlabel('Time (s)')
        ylabel('\Delta F / F')
        xticks(timeticks)
        xticklabels(xlabels)
        axis tight
        
    subplot(4,3,5)
        plot(zero_and_max(cns4s.intensity(1:pt_max,S4(i)).'),'color',cmap(i,:)), hold on
        u = fill(x,y,cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none'; hold on
        ylim([-0.2,max(zero_and_max(cns4s.intensity(1:pt_max,S4(i)).'))+0.1])
        spikes = cns4s.spikes(1:pt_max,S4(i)).'; 
        spktimes = find(spikes==1); 
        Nspk = length(spktimes);
        for ii = 1:Nspk
            xs = [spktimes(ii),spktimes(ii)];
            ys = [-0.15,-0.05];
            plot(xs,ys,'color','k'), hold on
        end
        xlabel('Time (s)')
        xticks(timeticks)
        xticklabels(xlabels)
        axis tight
        
    subplot(4,3,6)
        plot(zero_and_max(cns5s.intensity(1:pt_max,S5(i)).'), 'color',cmap(i,:)), hold on 
        u = fill(x,y,cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        ylim([-0.2,max(zero_and_max(cns5s.intensity(1:pt_max,S5(i)).'))+0.1])
        spikes = cns5s.spikes(1:pt_max,S5(i)).'; 
        spktimes = find(spikes==1); 
        Nspk = length(spktimes);
        for ii = 1:Nspk
            xs = [spktimes(ii),spktimes(ii)];
            ys = [-0.15,-0.05];
            plot(xs,ys,'color','k'), hold on
        end
        xlabel('Time (s)')
        xticks(timeticks)
        xticklabels(xlabels)
        axis tight
        
    tstrt = 2000;
    tstp = 2500; 
    timeticks = linspace(1,tstp-tstrt+1,5);
    timeticks2 = linspace(tstrt,tstp,5);
    xlabels = cell(5,1);
    for k = 1:5
        xlabels(k) = num2cell(round(timeticks2(k)/30));
    end
    subplot(4,3,7)
        data = zero_and_max(cns2s.intensity.').'; 
        plot(data(tstrt:tstp,S2(i)), 'color',cmap(i,:)), axis tight,  hold on
        u = fill([1 1 500 500],[-0.2, max(data(tstrt:tstp,S2(i)))+0.1, max(data(tstrt:tstp,S2(i)))+0.1, -0.2],cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        spikes = cns2s.spikes(tstrt:tstp,S2(i)).'; 
        spktimes = find(spikes==1); 
        Nspk = length(spktimes);
        for ii = 1:Nspk
            xs = [spktimes(ii),spktimes(ii)];
            ys = [-0.15,-0.05];
            plot(xs,ys,'color','k'), hold on
        end
        xticks(timeticks)
        xticklabels(xlabels)
        xlabel('Time (s)')
        ylabel('\Delta F / F')
        
    subplot(4,3,8)
        data = zero_and_max(cns4s.intensity.').';
        plot(data(tstrt:tstp,S4(i)), 'color',cmap(i,:)), axis tight,  hold on
        u = fill([1 1 500 500],[-0.2, max(data(tstrt:tstp,S4(i)))+0.1, max(data(tstrt:tstp,S4(i)))+0.1, -0.2],cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        spikes = cns4s.spikes(tstrt:tstp,S4(i)).'; 
        spktimes = find(spikes==1); 
        Nspk = length(spktimes);
        for ii = 1:Nspk
            xs = [spktimes(ii),spktimes(ii)];
            ys = [-0.15,-0.05];
            plot(xs,ys,'color','k'), hold on
        end
        xticks(timeticks)
        xticklabels(xlabels)
        xlabel('Time (s)')
    subplot(4,3,9)
        data = zero_and_max(cns5s.intensity.').';
        plot(data(tstrt:tstp,S5(i)), 'color',cmap(i,:)), axis tight, hold on
        u = fill([1 1 500 500],[-0.2, max(data(tstrt:tstp,S5(i)))+0.1, max(data(tstrt:tstp,S5(i)))+0.1, -0.2],cmap(i,:));
        u.FaceAlpha = 0.15;
        u.EdgeColor = 'none';
        spikes = cns5s.spikes(tstrt:tstp,S5(i)).'; 
        spktimes = find(spikes==1); 
        Nspk = length(spktimes);
        for ii = 1:Nspk
            xs = [spktimes(ii),spktimes(ii)];
            ys = [-0.15,-0.05];
            plot(xs,ys,'color','k'), hold on
        end
        xticks(timeticks)
        xticklabels(xlabels)
        xlabel('Time (s)')
    
        bw = 0.3;
    subplot(4,3,10)
        %ISIs = simple_ISIs(cns2s.mspikes(i,:));
        %h = histogram(ISIs,'BinWidth',0.5,'Normalization','probability'); hold on
        %h.FaceColor = cmap(i,:);
        %h.FaceAlpha = 0.8;
        ISIs2 = simple_ISIs(cns2s.spikes(:,S2(i)));
        g = histogram(ISIs2,'BinWidth',bw,'Normalization','probability');
        g.FaceColor = cmap(i,:);
        g.FaceAlpha = 0.3; 
        [N,~] = histcounts(ISIs2,'BinWidth',bw,'Normalization','Probability');
        whbin = find(N == max(N));
        modefreq = whbin * bw;
        title(['Mode = ',num2str(modefreq),'Hz'])
        xlim([0,10])
        xlabel('Hz bin')
        ylabel('Probability')
    subplot(4,3,11)
        %ISIs = simple_ISIs(cns4s.mspikes(i,:));
        %h = histogram(ISIs,'BinWidth',0.5,'Normalization','probability'); hold on
        %h.FaceColor = cmap(i,:);
        %h.FaceAlpha = 0.8;
        ISIs2 = simple_ISIs(cns4s.spikes(:,S4(i)));
        g = histogram(ISIs2,'BinWidth',bw,'Normalization','probability');
        g.FaceColor = cmap(i,:);
        g.FaceAlpha = 0.3;
        [N,~] = histcounts(ISIs2,'BinWidth',bw,'Normalization','Probability');
        whbin = find(N == max(N));
        modefreq = whbin * bw;
        title(['Mode = ',num2str(modefreq),'Hz'])
        xlim([0,10])
        xlabel('Hz bin')
    subplot(4,3,12)
        %ISIs = simple_ISIs(cns5s.mspikes(i,:));
        %h = histogram(ISIs,'BinWidth',0.5,'Normalization','probability'); hold on
        %h.FaceColor = cmap(i,:);
        %h.FaceAlpha = 0.8; hold on
        ISIs2 = simple_ISIs(cns5s.spikes(:,S5(i)));
        g = histogram(ISIs2,'BinWidth',bw,'Normalization','probability');
        g.FaceColor = cmap(i,:);
        g.FaceAlpha = 0.3; 
        [N,~] = histcounts(ISIs2,'BinWidth',bw,'Normalization','Probability');
        whbin = find(N == max(N));
        modefreq = whbin * bw;
        title(['Mode = ',num2str(modefreq),'Hz'])
        xlim([0,10])
        xlabel('Hz bin')
end


%% Session-wide average ISIs

cmap = jet(8);
bw = 0.3;
figure, hold on

subplot(2,3,1)
purkinje_artscape(cns2s.mask)
title('Session 2')

subplot(2,3,2)
purkinje_artscape(cns4s.mask)
title('Session 4')
subplot(2,3,3)
purkinje_artscape(cns5s.mask)
title('Session 5')

subplot(2,3,4)
ISIs = [];
for k = 1:cns2s.n_cells
    ISI = simple_ISIs(cns2s.spikes(:,k));
    ISIs = categ(1,ISI,ISIs);
end
h = histogram(ISIs,'BinWidth',bw,'Normalization','Probability');
[N,~] = histcounts(ISIs,'BinWidth',bw,'Normalization','Probability');
whbin = find(N == max(N));
modefreq = whbin * bw;
title(['Mode = ',num2str(modefreq),'Hz'])
h.FaceColor = cmap(1,:);
h.FaceAlpha = 0.3;
xlim([0,10])
xlabel('Frequency (Hz)')
ylabel('Probability')

subplot(2,3,5)
ISIs = [];
for k = 1:cns4s.n_cells
    ISI = simple_ISIs(cns4s.spikes(:,k));
    ISIs = categ(1,ISI,ISIs);
end
h = histogram(ISIs,'BinWidth',bw,'Normalization','Probability');
[N,~] = histcounts(ISIs,'BinWidth',bw,'Normalization','Probability');
whbin = find(N == max(N));
modefreq = whbin * bw;
title(['Mode = ',num2str(modefreq),'Hz'])
h.FaceColor = cmap(2,:);
h.FaceAlpha = 0.25;
xlim([0,10])
xlabel('Frequency (Hz)')

subplot(2,3,6)
ISIs = [];
for k = 1:cns5s.n_cells
    ISI = simple_ISIs(cns5s.spikes(:,k));
    ISIs = categ(1,ISI,ISIs);
end
h = histogram(ISIs,'BinWidth',bw,'Normalization','Probability');
[N,~] = histcounts(ISIs,'BinWidth',bw,'Normalization','Probability');
whbin = find(N == max(N));
modefreq = whbin * bw;
title(['Mode = ',num2str(modefreq),'Hz'])
h.FaceColor = cmap(3,:);
h.FaceAlpha = 0.2;
xlim([0,10])
xlabel('Frequency (Hz)')



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
m_cas5 = mean(cns5s.intensity.',1);
m_cas5 = zero_and_max(m_cas5);

m_cas2 = mean(cns2s.intensity.',1);
m_cas2 = zero_and_max(m_cas2);

m_cas4 = mean(cns4s.intensity.',1);
m_cas4 = zero_and_max(m_cas4); 


Nkk = 5;
inputpath = '/Volumes/carey/leonard.dupont/TM IMAGING FILES/voluntary locomotion/MC318/S2';
speedkk2 = imbased_konkat(im2,tm2,inputpath,Nkk);
tm_speedMs2  = speedkk2.trials_15_16_1_27_29;

Nkk = 5;
inputpath = '/Volumes/carey/hugo.marques/LocomotionExperiments/RT Self-paced/Imaging/TM IMAGING FILES/voluntary locomotion/MC318/S4';
speedkk4 = imbased_konkat(im4,tm4,inputpath,Nkk);
tm_speedMs4  = speedkk4.trials_1_21_24_25_31;

Nkk = 5;
inputpath = '/Volumes/carey/hugo.marques/LocomotionExperiments/RT Self-paced/Imaging/TM IMAGING FILES/voluntary locomotion/MC318/S5';
speedkk5 = imbased_konkat(im5,tm5,inputpath,Nkk);
tm_speedMs5  = speedkk5.trials_1_25_27_51;

clear speedkk2
clear speedkk4
clear speedkk5

%% Producing gaussian-filtered speed traces 

%removing negative speed values
tm_speedMs2(tm_speedMs2 < 0) = 0;
tm_speedMs4(tm_speedMs4 < 0) = 0;
tm_speedMs5(tm_speedMs5 < 0) = 0;

%gaussian filtering
tm_speedMs2 = smoothdata(tm_speedMs2,'gaussian',5000);
figure, plot(tm_speedMs2)
tm_speedMs4 = smoothdata(tm_speedMs4,'gaussian',5000);
figure, plot(tm_speedMs4)
tm_speedMs5 = smoothdata(tm_speedMs5,'gaussian',5000); %smoothed version
figure, plot(tm_speedMs5)

%%
% To define speed categories, given that it is the same mouse each time, we
% are going to use a concatenated speed trace.
ncat = 3; 
l_ALL = length(tm_speedMs2)+length(tm_speedMs4)+length(tm_speedMs5);
tm_speedALL = zeros(l_ALL,1);
tm_speedALL(1:length(tm_speedMs2),1) = tm_speedMs2;
tm_speedALL((length(tm_speedMs2)+1):(length(tm_speedMs2)+length(tm_speedMs4)),1) = tm_speedMs4;
tm_speedALL((length(tm_speedMs2)+length(tm_speedMs4)+1):l_ALL,1) = tm_speedMs5;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
speedcats = define_speedcats(tm_speedALL,ncat);
% - - - - - - - - - - - - OR (individual speedcats) - - - - - - - - - - - - 
speedcats2 = define_speedcats(tm_speedMs2,ncat);
speedcats4 = define_speedcats(tm_speedMs4,ncat);
speedcats5 = define_speedcats(tm_speedMs5,ncat);
% - - - - - - - - - Mind the input speedcats - - - - - - - - - - - - - - -
speedlabels2 = assign_speedlabels(tm_speedMs2,speedcats2);
speedlabels4 = assign_speedlabels(tm_speedMs4,speedcats4);
speedlabels5 = assign_speedlabels(tm_speedMs5,speedcats5);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
xc2 = get_speedregime_boundaries(speedlabels2);
xc4 = get_speedregime_boundaries(speedlabels4);
xc5 = get_speedregime_boundaries(speedlabels5);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot_all_traces_w_speed(cns2s.intensity,xc2,speedlabels2,tm_speedMs2)
plot_all_traces_w_speed(cns4s.intensity,xc4,speedlabels4,tm_speedMs4)
plot_all_traces_w_speed(cns5s.intensity,xc5,speedlabels5,tm_speedMs5)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
%% Acceleration, decceleration and constant speed
% Here we define three acceleration regimes and see how these relate to
% cell activity. 

accs4 = get_acceleration_from_speed(tm_speedMs4,tm4.time(1:length(tm_speedMs4)));
mm = mean(accs4);
accs4 = accs4 - mm;  %center
ups = max(accs4);
downs = min(accs4); 
sigmaup = 0.02;
sigmadown = -0.02;

figure, plot(accs4,'color','k'), axis tight, hold on
hline(sigmaup), hold on
hline(sigmadown), 
xlabel('Time')
ylabel('Acceleration')
hold off 


% Now we classify data 
acclabels = zeros(length(accs4),1);
for k = 1:length(accs4)
    pt = accs4(k);
    if (pt<sigmaup) && (pt>sigmadown)
        acclabels(k) = 2; %small fluctuations
    elseif pt > sigmaup
        acclabels(k) = 3; %noteworthy acceleration
    elseif pt < sigmadown
        acclabels(k) = 1; %noteworthy decceleration
    end
end
%%
acccolors = zeros(length(accs4),3);
for k = 1:length(acclabels)
    switch acclabels(k)
        case 1
            acccolors(k,:) = [0 0.3 0.9];
        case 3
            acccolors(k,:) = [1 0.4 0.3];
        case 2
            acccolors(k,:) = [0.8 0.8 0.8];
    end
end
xcacc = get_speedregime_boundaries(acclabels);
fluo_x_accelerationpeaks(cns4.intensity,accs4,acccolors,xcacc)
%% PMI for mean calcium trace

m_cas4 = mean(cns4.intensity,2).'; 

naccat = 3;
colormap = hot(naccat);

nfluocats = 20;
fluocats = define_fluocats(m_cas4.',nfluocats);
fluolabels = assign_fluolabels(m_cas4.',fluocats);

Tim = length(m_cas4);
Ttm = length(accs4);
nsynchbins = 10;
mergeF = Ttm/Tim; 

binboundaries = zeros(Tim,1); 
for k = 1:Tim
    binboundaries(k) = floor((k-1)*mergeF + 1);
end
aside = Ttm - binboundaries(end);
if aside ~=0
    binboundaries(end+1) = binboundaries(end)+aside;
end

Nimbins = max(fluolabels);
Ntmbins = max(acclabels);
MImat = zeros(Nimbins,Ntmbins);

for k = 1:Tim
    imlbl = fluolabels(k);
    start = binboundaries(k);
    stop = binboundaries(k+1);
    acclbl = round(mean(acclabels(start:stop)));
    MImat(imlbl,acclbl) = MImat(imlbl,acclbl) + 1;
end

MImat = MImat / Tim; 

for imlbl = 1:Nimbins
    Pim = length(find(fluolabels == imlbl)) / length(fluolabels);
    for acclbl = 1:Ntmbins
        Ptm = length(find(acclabels == acclbl)) / length(acclabels);
        MImat(imlbl,acclbl) = log(MImat(imlbl,acclbl)/(Pim*Ptm));
    end
end

figure, imagesc(MImat)
colorbar
ylabel('Fluorescence bins')
xlabel('Acceleration regime')
xticks([1,2,3])
xticklabels({'Decceleration','Ground state','Acceleration'})


%%
nfluocats = 20;
Nimbins = nfluocats;
clear PMI
for roi = 1:cns4.n_cells
    
    MImat = zeros(Nimbins,Ntmbins);
    fluocats = define_fluocats(cns4.intensity(:,roi),nfluocats);
    fluolabels = assign_fluolabels(cns4.intensity(:,roi),fluocats);

    for k = 1:Tim
        imlbl = fluolabels(k);
        start = binboundaries(k);
        stop = binboundaries(k+1);
        acclbl = round(mean(acclabels(start:stop)));
        MImat(imlbl,acclbl) = MImat(imlbl,acclbl) + 1;
    end

    MImat = MImat / Tim; 
    
    for imlbl = 1:Nimbins
        Pim = length(find(fluolabels == imlbl)) / length(fluolabels);
        for acclbl = 1:Ntmbins
            Ptm = length(find(acclabels == acclbl)) / length(acclabels);
            if Pim>0 && Ptm>0
                MImat(imlbl,acclbl) = MImat(imlbl,acclbl)/(Pim*Ptm);
            else
                MImat(imlbl,acclbl) = 0;
            end
        end
    end
    
    mMI = max(MImat(:));
    MImat = MImat / mMI;
    colormap('Parula')
    subplot(10,19,roi), imagesc(MImat), axis off, title(num2str(roi))
    PMI.(['roi_',num2str(roi)]) = MImat;
    
end
suptitle('Cell-resoluted PMI between acceleration and fluorescence')

%%
PMIcoord = zeros(cns4.n_cells,3);
start = 17 ;
stop = nfluocats;
for roi = 1:cns4.n_cells
    MImat = PMI.(['roi_',num2str(roi)]);
    for j = 1:3
       PMIcoord(roi,j) = mean(MImat(start:stop,j))/mean(MImat(:));
    end
end

subplot(1,2,1)
scatter3(PMIcoord(:,1),PMIcoord(:,2),PMIcoord(:,3),'filled'), box off
xlabel('Decceleration')
ylabel('Ground state')
zlabel('Acceleration')

K = 3;
indx = kmeans(PMIcoord,K);
cmap = parula(K);
subplot(1,2,2), hold on
for roi = 1:cns4.n_cells
    ii = indx(roi);
    scatter3(PMIcoord(roi,1),PMIcoord(roi,2),PMIcoord(roi,3),[],cmap(ii,:),'filled'), grid on
    text(PMIcoord(roi,1),PMIcoord(roi,2),PMIcoord(roi,3),num2str(roi),'color',cmap(ii,:))
end

xlabel('Decceleration')
ylabel('Ground state')
zlabel('Acceleration')
hold off 
%%
figure, hold on
for k = 1:K
    cl = find(indx == k);
    L = length(cl);
    mintensity = zeros(length(cns4.intensity(:,1)),L);
    mPMI = zeros(20,3);
    for i = 1:L
        mPMI = mPMI + PMI.(['roi_',num2str(cl(i))]);
        mintensity(:,i) = zero_and_max(cns4.intensity(:,cl(i)).').';  
    end
    mPMI = mPMI/L;
    subplot(1,K,k), colormap('Gray'), imagesc(mPMI), box off 
    colorbar 
    if k == 1
        ylabel('Fluorescence bins')
        xlabel('Acceleration regime')
        xticks([1,2,3])
        xticklabels({'Decceleration','Ground state','Acceleration'})
    else
        xticks([])
        yticks([])
    end
    title(['Cluster ',num2str(k)])
    clintensity.(['cluster_',num2str(k)]) = mintensity; 
end
hold off 

%%
figure, hold on
offset = 0;
cvmat = zeros(K,K);
for k = 1:K
    inte = clintensity.(['cluster_',num2str(k)]);
    inte = mean(inte,2);
    inte = inte - min(inte); 
    clintensity.(['cluster_',num2str(k)]) = inte; 
    plot(inte,'color',cmap(k,:),'LineWidth',1)
end

%%
inte1 = clintensity.(['cluster_',num2str(1)]);
inte2 = clintensity.(['cluster_',num2str(2)]);
inte3 = clintensity.(['cluster_',num2str(3)]);
sortedvals = zeros(length(inte1),1); 
 
figure, hold on
time = linspace(1,length(accs4),length(inte1)); 
pcateg = 1 ;
start = 1;
stop = 1;
down = 0;
for j = 1:K
    inte = clintensity.(['cluster_',num2str(j)]);
    for k = 1:length(inte1)
        m = inte(k);
        categ = 1;
        while m > fluocats(categ)
            categ = categ + 1;
        end
        if categ ~= pcateg 
            xl = [start, start, stop, stop];
            yl = [down, -0.05 + down, -0.05+down, down];
            b = fill(xl*100,yl,cmap(j,:)*categ/20); hold on,
            b.FaceAlpha = 1;
            b.EdgeColor = 'none';
            start = k;
            stop = k + 1;
            categ = maxi; 
        else
            stop = stop + 1;
        end

    end
    down = down - 0.05;
end
plot(time,inte1,'color',cmap(1,:),'LineWidth',1), hold on
plot(time,inte2,'color',cmap(2,:),'LineWidth',1), hold on
plot(time,inte3,'color',cmap(3,:),'LineWidth',1), hold on

for k = 1:length(xcacc)-1
    xplus = xcacc(k+1);
    xminus = xcacc(k);
    color = acccolors(xplus-1,:);
    
    xl = [xminus, xminus, xplus, xplus];
    yl = [-0 1.1 1.1 -0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.2;
    a.EdgeColor = 'none';
end



ylabel('\Delta F/F')

axis tight 

%%

conditioned_fluovals = zeros(K,naccat);
dsacclabels = zeros(Tim,1); 
for j = 1:K
    cl = find(indx == j);
    L = length(cl);
    nacclbl = zeros(3,1);
    acc = [];
    cst = [];
    decc = [];

    for i = 1:L
        theintensity = zero_and_max(cns4.intensity(:,cl(i)).'); 
        for k = 1:Tim
            fluoval = theintensity(k); 
            start = binboundaries(k);
            stop = binboundaries(k+1);
            acclbl = round(mean(acclabels(start:stop)));
            dsacclabels(k) = acclbl;
            if fluoval > 0.3
                conditioned_fluovals(j,acclbl) = conditioned_fluovals(j,acclbl) + fluoval; 
                nacclbl(acclbl) =  nacclbl(acclbl) + 1;
                
            end   
            switch acclbl
                    case 1
                        decc(end+1) = fluoval;
                    case 2
                        cst(end+1) = fluoval;
                    case 3
                        acc(end+1) = fluoval; 
            end 
 
        end
    end
    
    figure, title(['Cluster ',num2str(j)]), hold on
    h1 = histogram(decc,'binwidth',0.005,'Normalization','countdensity');
    h1.FaceColor = [0 0.3 0.9];
    h1.FaceAlpha = 0.5;
    
    h2 = histogram(cst,'binwidth',0.005,'Normalization','countdensity');
    h2.FaceColor =  [0.8 0.8 0.8];
    h2.FaceAlpha = 0.5;
    
    h3 = histogram(acc,'binwidth',0.005,'Normalization','countdensity');
    h3.FaceColor = [1 0.4 0.3];
    h3.FaceAlpha = 0.5;
    
    hold off 
   
    for k = 1:naccat
        conditioned_fluovals(j,k) = conditioned_fluovals(j,k)/nacclbl(k);
    end

    
end

figure, subplot(1,3,1), bar(conditioned_fluovals(1,:))
subplot(1,3,2), bar(conditioned_fluovals(2,:))
subplot(1,3,3), bar(conditioned_fluovals(3,:))


%%
figure, hold on
cmap = parula(K);
S = size(cns4.mask{1,1});
for k = 1:K
    cl = find(indx == k);
    L = length(cl);
    c = cmap(k,:);
    full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
    for i = 1:L
        roi = cl(i);
        posi = cns4.centroid{1,roi};
        x = posi(1);
        y = posi(2);
        text(x,y,num2str(roi),'Color',c)
        I = cns4.mask{1,roi};
        h = imshow(full); hold on 
        set(h, 'AlphaData', I*0.35) , hold on
    end
end
%%

clear opts
opts.Nmax = 10;
opts.Nmin = 2;
opts.epsilon = 10e-4;
opts.framepath = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/registration_templateS4.tif';


clear grphcs
grphcs.dendrogram = 1;
grphcs.orgdistMAT = 1;
grphcs.orgchanceMAT = 1;
grphcs.porgchanceMAT = 1;
grphcs.clusteredlandscape = 1;

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
synchresults = synchrony_clustering(cns4s,opts,grphcs); 



%% END OF ACCELERATION x FLUO PMI ANALYSIS
%% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
%% Now let's move onto some actual fluo x speed analysis
%     Synchrony
%     Using the acceleration
%     Fluorescence intensity

tic
wsize = 6; %points
dw = 1;


% this computation is going to be huge. 
N = cns4s.n_cells;
Nw = floor(length(cns4s.spikes(:,1))/dw);
aside = rem(length(cns4s.spikes(:,1)),dw);
ssynchrony = zeros(Nw,1);

for k = 1:Nw-1
    start = (k-1)*dw + 1;
    if start + wsize - 1 < length(cns4s.spikes(:,1))
        stop = start + wsize - 1 + (k==Nw)*aside;
    else
        disp('Finished')
        break
    end
    Qt = zeros(N*(N-1)/2,1);
    c = 1;
    for ii = 1:N
        ispk = cns4s.spikes(start:stop,ii);
        for jj = 1:N
            if jj > ii
            jspk = cns4s.spikes(start:stop,jj);
            [Q,~,~,~] = est_trace_synchrony(ispk,jspk);
            Qt(c) = Q;
            c = c + 1;
            end
        end
    end
    ssynchrony(k) = mean(Qt);
    if rem(k,100) == 0 
        disp(['Window ',num2str(k),' out of ',num2str(N),'.'])
    end
end
clear c

%%
sm_synchrony = smoothdata(ssynchrony,'gaussian',10); 
figure, plot(sm_synchrony)

%% PSTH
% Mystery with this plot is time alignement of subplots
% 1 : rasterplot with acceleration colorcode
% 2 : psth-like plot
% 3 : corresponding synchrony = f(t) (should look alike)
ah1 = subplot(7,1,1:5);
height = 1;
spacing = 0.5;
N = cns4s.n_cells;
im_max = 145*30;
tm_max = 145*3000; 
xcacc(end) = tm_max;

mH = roi*height + roi*spacing;
for k = 1:length(xcacc)-1
    xplus = xcacc(k+1);
    xminus = xcacc(k);
    color = acccolors(xplus-1,:);
    
    xl = [xminus, xminus, xplus, xplus]*1/3000;
    yl = [0  mH mH 0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.6;
    a.EdgeColor = 'none';
end


for roi = 1:N
        
   estimated = cns4s.spikes(1:im_max,roi);
   spktimes = find(estimated==1); %each round, we get the spike times 
   y = [(roi-1)*height + roi*spacing , roi*height + roi*spacing]; %we prepare 2 coordinates in y (we're gonna draw a vertical line: spike)
        
   for spk = 1:length(spktimes) %for each spike from this roi, we draw vertical lines at the right x coordinates
       x = [spktimes(spk), spktimes(spk)]*1/30.0003;
       plot(x,y,'color','k'), hold on
   end
        
end

set(gca,'TickLength',[0.001,0])
box off 
xlim([1,145])
ylim([0,mH])
ylabel('Units')
xticks([])
yticks([])

ah2 = subplot(7,1,6);
sumSPK = zeros(length(cns4s.intensity(1:im_max,1)),1);
for k = 1:length(cns4s.spikes(1:im_max,1))
    spkthistime = cns4s.spikes(k,:);
    sumSPK(k) = sum(spkthistime);
end

for k = 1:length(xcacc)-1
    xplus = xcacc(k+1);
    xminus = xcacc(k);
    color = acccolors(xplus-1,:);
    
    xl = [xminus, xminus, xplus, xplus]*1/3000;
    yl = [0 max(sumSPK) max(sumSPK) 0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.6;
    a.EdgeColor = 'none';
end

xb = linspace(1,length(sumSPK)/30,length(sumSPK));
clear spkthistime
b = bar(xb,sumSPK); hold on 
b.FaceColor = 'k';
%b.BarWidth = 1;
b.FaceAlpha = 1;
b.EdgeColor = 'k';
b.EdgeAlpha = 1;

xlim([1,145])
xticks([])
ylabel('N_{spikes}')

ah3 = subplot(7,1,7);
for k = 1:length(xcacc)-1
    xplus = xcacc(k+1);
    xminus = xcacc(k);
    color = acccolors(xplus-1,:);
    
    xl = [xminus, xminus, xplus, xplus]*1/3000;
    yl = [0 1 1 0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.6;
    a.EdgeColor = 'none';
end

x = linspace(1,im_max/30,im_max);
plot(x,sm_synchrony(1:im_max),'color','k'), hold on


tt = linspace(1,im_max/30,10);

lbls = cell(length(tt),1);
for j = 1:length(tt)
    lbls(j) = num2cell(round(tt(j)));
end
ylabel('Synchrony')
xticks(tt)
xticklabels(lbls)
xlabel('Time (s)')
xlim([1,145])
box off
hold off 

linkaxes([ah1,ah2,ah3],'x')

%% EVENT AMPLITUDE AND CELL SYNCHRONY
% It seems that because of scattering in the 1P microscope, events that
% correspond to huge synchronisations are reflected by high-fluorescence
% events in the calcium traces. Here we quantify this phenomenon.

for i = 1:cns4.n_cells
    spikes = cns4.spikes(:,i);
    calcium = (cns4.intensity(:,i) - min(cns4.intensity(:,i)))/max(cns4.intensity(:,i)); %assuming flat baseline
    spktimes = find(spikes == 1);
    eventheight = zeros(length(cns4.spikes(:,i)),1);
    for k = 1:length(spktimes)
        t = spktimes(k);
        eventheight(t) = mean(calcium(t:t+2));
    end
    cns4.eventheight(:,i) = eventheight; 
end



clear spkamp_coord
c = 1;
for k = 1:length(cns4.intensity)
    spkunits = find(cns4.spikes(k,:) == 1);
    Nunits = length(spkunits);
    if Nunits > 1
        posvect = zeros(Nunits,2);
        for i = 1:Nunits
            roi = spkunits(i);
            posvect(i,:) = cns4.centroid{1,roi};
        end
        xc = mean(posvect(:,1));
        yc = mean(posvect(:,2));
        maxdist = 0;
        for i = 1:Nunits
            x = posvect(i,1);
            y = posvect(i,2);
            mdist = sqrt((x-xc)^2 + (y-yc)^2); 
            if mdist > maxdist
                maxdist = mdist;
            end
        end
 
    else
        maxdist = 0;
    end
    spkamp_coord(1,c) = mean(cns4.eventheight(k,:));
    spkamp_coord(2,c) = Nunits;
    spkamp_coord(3,c) = maxdist; 
    c = c + 1;
end
%%
varnames = {'Event Amplitude','Number of synchronous units','Max distance to centroid'};
c = 1;
for k = 1:3
    a = spkamp_coord(k,:);
    for j = 1:3
        b = spkamp_coord(j,:);
        if k~=j
            subplot(2,3,c), hold on
            s = scatter(a,b,3,cmap(k,:),'filled'); hold on
            s.MarkerFaceAlpha = 1;
            s.MarkerEdgeColor = cmap(k,:);
            s.MarkerEdgeAlpha = 0.5;
            xlabel(varnames{k})
            ylabel(varnames{j})
            c = c+1;
        end
    end
end

subplot(3,3,7)
            
x = spkamp_coord(2,:); %Nunits sharing the spike  
y = spkamp_coord(1,:);
z = spkamp_coord(3,:); 

zz = find(x==0);
x(zz(1:end-1)) = [];
y(zz(1:end-1)) = [];
z(zz(1:end-1)) = [];

s = scatter3(x,y,z,3,cmap(1,:),'filled'); 
s.MarkerFaceAlpha = 1;
s.MarkerEdgeColor = cmap(1,:);
s.MarkerEdgeAlpha = 0.5;
xlabel('Number of synchronous cells')
ylabel('Event amplitude (\Delta F/ F)')
zlabel('Maximum distance between synchronous cells')
%%
subplot(1,2,1)
x = spkamp_coord(2,:); %Nunits sharing the spike  
y = spkamp_coord(1,:); %amplitudes
zz = find(x==0);
x(zz(1:end-1)) = [];
y(zz(1:end-1)) = [];
s = scatter(x,y,3,[0.767816699258443   0.516358225453477   0.506171068243531],'filled'); hold on
s.MarkerFaceAlpha = 1;
s.MarkerEdgeColor = [0.767816699258443   0.516358225453477   0.506171068243531];
s.MarkerEdgeAlpha = 0.5;
xlabel('Number of synchronous cells')
ylabel('Event amplitude (\Delta F/ F)')


X = cat(1,ones(1,length(x)),x);
b = X'\y';
y_prd= X'*b; %predicted y based on regression 
plot(x,y_prd,'--','color','k'), hold on

Rsquare = 1 - sum((y - y_prd.').^2)/sum((y-mean(y)).^2);
title(['r^{2} = ',num2str(Rsquare)])
ylim([0,max(y)])
box off 


% - - - - - - - amplitude categories - - - - - - - -
nampcats = 50;
maxy = max(y);
diffcat = maxy/nampcats ; %mins is 0, fluorescence rarely is 0
ampcats = zeros(nampcats,1);
for i = 1:nampcats
   ampcats(i) = i* diffcat;
end

amplabels = zeros(1,length(y));
for k = 1:length(y)
    ampl = y(k);
    categ = 1;
    while ampl > ampcats(categ)
        categ = categ + 1;
    end
    amplabels(k) = categ;
end

% - - - - - - - same for Nunits - - - - - - -
nucats = 50;
maxx = max(x);
diffcat = maxx/nucats;
ucats = zeros(nucats,1);
for i = 1:nucats
   ucats(i) = i*diffcat;
end

ulabels = zeros(1,length(x));
for k = 1:length(x)
    nunits = x(k);
    categ = 1;
    while nunits > ucats(categ)
        categ = categ + 1;
    end
    ulabels(k) = categ;
end

% - - - - - - - Now PMI - - - - - - 
MImat = zeros(nampcats,nucats);
for k = 1:length(x)
    albl = amplabels(k);
    ulbl = ulabels(k);
    MImat(albl,ulbl) =  MImat(albl,ulbl) + 1;
end

MImat = MImat / length(x);

for aa = 1:nampcats
    Paa = length(find(amplabels == aa)) / length(amplabels);
    for uu = 1:nucats
        Puu = length(find(ulabels == uu)) / length(ulabels);
        MImat(aa,uu) =  log(MImat(aa,uu) /(Paa*Puu));
    end
end

subplot(1,2,2),colormap('Pink'), imagesc(MImat)
title('Point-Mutual information')
box off

xlabel('Synchronous-unit bin')
ylabel('Amplitude bin')

suptitle('Event amplitude as a function of units sharing the spike')

%%
cmap = pink(25); 
subplot(1,2,1)
h1 = histogram(x,'binwidth',1,'Normalization','Probability');
%h1.FaceColor = [0.980, 0.501, 0.501];
h1.FaceColor = cmap(8,:);
h1.FaceAlpha = 0.6;
xlabel('Number of synchronous cells')
ylabel('Probability')
box off 

subplot(1,2,2)
h2 = histogram(y,'binwidth',0.005,'Normalization','Probability');
xlim([1e-4,0.1])
%h2.FaceColor = [0.980, 0.694, 0.501];
h2.FaceColor = cmap(20,:);
h2.FaceAlpha = 0.6;
xlabel('Event amplitude (\Delta F /F)')
ylabel('Probability')
box off 

suptitle('Probability distributions')
%%

clear cl
cl.K = 5;
cl.runs = 2500 ;
cl.spatialplot = 1;
cl.frame_path = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/registration_template.tif';
cl.normalise = 0;
cl.p = 0.8;
cl.Z = 5;
cl.Nmin = 3;
cl.Nmax = 10; 

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
activity_clustersK4 = hybrid_clusteringDB(cns4,cl);
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

%%

clear opts4
opts4.framepath = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/registration_templateS4.tif' ;
opts4.wsize = 1;
opts4.Nmax = 10;
opts4.Nmin = 1;
opts4.min_units = 5;

clear grphcs
grphcs.dendrogram = 1;
grphcs.orgdistMAT = 1;
grphcs.orgchanceMAT = 1;
grphcs.porgchanceMAT = 1;
grphcs.clusteredlandscape = 1;

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
synchresultS4 = synchrony_clustering(cns4,opts4,grphcs); 
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

%%



