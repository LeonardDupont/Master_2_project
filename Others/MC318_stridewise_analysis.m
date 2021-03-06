%% MC318 - S4, stridewise analysis of calcium-imaging data
% leonard.dupont@ens.fr
% 
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
%           and align it to the right time.
%
% The extraction pipeline was devised by Diogo Duarte. 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%   We will first do this for the mean fluorescence trace and then go for a
%   cell-resoluted analysis. We also extract spikes.  
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%       TABLE OF CONTENTS                                      line
%
%   I   Whole session fluorescence extraction, deconvolution    37
%   II  Stride extraction                                       73
%   III Stride normalisation and storage                        116
%   IV  Bulk fluorescence extraction + normalisation            229
%   V   Cell-specific fluorescence extraction + norm.           297
%
%   VI  PLOTS                                                   388
%       a    Stride-locked bulk                                 396
%       b    Stride-locked cell matrices                        469
%       c    Stride-locked cell matrices (reorganised kmeans)   505
%       d    Stride-locked cluster averages (metakmeans)        590
%       e    Test of cluster significance (random splitting)    705
%       f    Evolution of variance throughout the stride        843
%       g    Stride-locked single-cell fluo.                    894
%       h    Amplitude of stride-locked modulation              959
%       i    Stride-locked complex-spike probabilities          1045
%

%% Colors, paths, others  - please compile before anything else!  

sessiontifffile = ''; 
tracking_dir = 'Z:\leonard.dupont\TM PROCESSED FILES\voluntary locomotion'; 
tdmsdata = 'Z:\leonard.dupont\TM SESSION FILES\voluntary locomotion\MC318\S4'; 
fstm = 300; %sampling frequency treadmill (Hz)
fsim = 30; %sampling frequency imaging (Hz)

%number of points to bin the stride and fluorescence
stdpts = 10; 
stdall = stdpts; 

% 1 : Front Right -- 2 : Hind Right -- 3 : Front Left -- 4 : Hind Left
paw_dictionary = {'FR','HR','FL','HL'}; 
paw_colors = [0.976, 0.223, 0.243 ; 0.972, 0.266, 0.831 ; 0.223, 0.388, 0.976 ; 0.223, 0.866, 0.976];



%% I - Extracting fluorescence from the whole session. 
% First, we build a cn-like struct with fluorescence data from the whole
% session. Fluorescence is extracted using the usual masks. 

wholeS4 = bigread2(sessiontiffile); %import tiff 
[~,~,t] = size(wholeS4); 

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

% . . . . . . . . deconvolution  (will be used after) . . . . . . . . . . . 
ops.fs = 30;
ops.recomputeKernel =0;
ops.sensorTau = 0.7;
ops.estimateNeuropil = 0;
ops.deconvType = 'L0'; 
%ops.deconvType = 'OASIS';
threshold = 1000;
[~, ~, ~, cnfull] = get_spikes_from_calcium_traces(cnfull.intensity, ops, threshold, cnfull, []);


%% II Extract strides

% loading imaging, tracking tdms as well as tracking files 
tracks = load(fullfile(tdmsdata,'tracking_data.mat')); 
tracks = tracks.tracks; 
im = load(fullfile(tdmsdata,'imaging_data.mat'));
im = im.im; 

% There we get all the .mat files that contain stride information 
filePattern = fullfile(tracking_dir,'*.mat');
Ntrials = length(dir(filePattern)); 
trials = dir(filePattern); 


% Unfortunately, these are not listed in the right order, so we just need
% to solve this using a secondary, indexed vector (called 'order' here). 
order = zeros(Ntrials,1); 
disorder = zeros(Ntrials,1); 
for k = 1:Ntrials
    name = trials(k).name; 
    underscores = find(name == '_'); 
    nb = name(underscores(end-1)+1:underscores(end)-1); 
    disorder(k) = str2double(nb);
end
for k = 1:Ntrials
    minidx = find(disorder == min(disorder));
    disorder(minidx) = 1e3; 
    order(minidx) = k; 
end

%% Separating trials and finding onset times based on tdms vectors
% This is both for imaging and tracking, although we won't need
% trial-based information for fluorescence since we're using a file with
% all concatenated sessions. 

tmidx = find_trial_onset(tracks.time,fstm); 


%% III  Isolating strides in groups of 'stdall' points
% Here, we travel through all the session trials and extract stride points
% (x values) as well as corresponding timestamps in the session. Everything
% is stored in a struct 'isolated_strides' with subfields for different
% paws. For each paw, we then have a matrix 'isolated' with datapoints of
% every detected stride for the limb in the session and a matrix 
% 'timestamps' with corresponding time values. 


% To sample the step cycle, we use three reference points : stance offset (1),
% swing (2) and stance onset (3). We add 4 points between (1) and (2) and 3
% points between (2) and (3). This can be changed hereunder of course. 


countit = [0 0 0 0]; %this is to correctly build the struct (see below)
% we count if this the first time we treat the paw in question 

clear isolated_strides


% we now loop through the trials 
for k = 1:Ntrials
    
    % We use the 'order' vector to read trials chronologically 
    realidx = find(order == k); 
    name = trials(realidx).name; 
    disp(['[stride extraction] Processing ', name, '.'])
    dtMAT = load([tracking_dir,'\',name]);
    
    % we keep stride reference points ('strides') and paw trajectories
    % ('paws') 
    strides = dtMAT.strides; 
    paws = dtMAT.trial_data.paws_x;
    
    % We use the output indices from find_trial_onset.m 
    if k < Ntrials
        ttime = tracks.time(tmidx(k):tmidx(k+1)); 
    else %if it is the last trial, then we go till the end
        ttime = tracks.time(tmidx(k):end);
    end
    
    %for the given trial, we now loop through the 4 paws
    for pw = 1:4
        coord = strides{pw}; %reference points of all strides for the paw
        Nstr = length(coord); %how many strides are there? 
        
        %allocation of our matrices 
        isolated = zeros(Nstr,stdall); 
        timestamps = zeros(Nstr,stdall); 
        
        %for the given trial and paw, we loop through the tracked strides 
        for st = 1:Nstr

            refs = coord(st,:); %reference points for this stride
            
            % Hereunder, we prepare evenly spaced sample points to keep
            % x-trajectory values between the reference points of the
            % stride (see explanation above). 
            samplepoints = zeros(stdall,1);
            samplepoints(1) = refs(1); 
            samplepoints(end) = refs(3);

            %first bit
            dt = (refs(3)-refs(1))/(stdpts1-1); %refs(1) to refs(2)
            for dd = 1:stdpts
                samplepoints(1+dd) = round(dd*dt) + refs(1); 
            end 
            

            % We can now extract values and store timestamps for the stride
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            isolated(st,:) = paws(samplepoints,pw);
            timestamps(st,:) = ttime(samplepoints);
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        end
        
        % once all strides of the paw have been done, we add the data to
        % the paw subfield. 
        
        if countit(pw) 
        %if it's not the first time we sort data for the paw we concatenate
        %the matrix with the previously stored ones. 
            isolated = cat(1,isolated_strides.([paw_dictionary{pw},'_paw']).isolated,isolated); 
            timestamps = cat(1,isolated_strides.([paw_dictionary{pw},'_paw']).timestamps,timestamps); 
        end
        
        
        isolated_strides.([paw_dictionary{pw},'_paw']).isolated = isolated; 
        isolated_strides.([paw_dictionary{pw},'_paw']).timestamps = timestamps;  
        countit(pw) = 1; 
    end
end

%% Plot the mean stride for each paw
figure, hold on
for pw = 1:4
    subplot(2,2,pw), hold on
    data = isolated_strides.([paw_dictionary{pw},'_paw']).isolated; 
    [x,~] = size(data);
    for k = 1:x
        plot(data(k,:),'Color',[0.8 0.8 0.8]), hold on
    end
    plot(mean(data,1),'LineWidth',1.5,'Color','k')
    xlabel('Time (datapoints)')
    ylabel('x coord.')
    axis tight, box off 
    title(paw_dictionary{pw})
end

%% IV  Isolating fluorescence in groups of 'stdall' points 
% Now that we have extracted strides for each trial in the session, we'd
% like to snip the corresponding fluorescence values, either from bulk
% fluorescence or for individual cells. Identically, we will sample
% (stdall) points of fluorescence using the stride timestamps and store
% everything in a structure, 'isolated_fluorescence' with paw-specific
% subfields and 'isolated' + 'timestamps' further subsubfields. 

clear isolated_fluorescence
mintensity = cnf.mintensity;

% just like for the previous struct 

% We will loop through paws. 
for pw = 1:4
  
  % First we extract the corresponding tracking data 
  isolated = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
  timestamps = isolated_strides.([paw_dictionary{pw},'_paw']).timestamps;
  
  %How many strides will we snip fluorescence for 
  [Nstr,~] = size(isolated);

  %preallocation
  isofluo = zeros(Nstr,stdall); 
  imtimestamps = zeros(Nstr,stdall); 
  
  % We now loop through each stride of this paw 
  for ss = 1:Nstr
      onset = timestamps(ss,1); %refs(1)
      offset = timestamps(ss,end); %refs(3) 
      
      % In this paragraph, we find the closest timestep in im.time for both
      % the stride onset and offset.
      whereon = abs(im.time - onset);
      imon = find(whereon == min(whereon)) - 1; 
      imon = imon(1); 
      whereoff = abs(im.time - offset);
      imoff = find(whereoff == min(whereoff)) + 1; 
      imoff = imoff(1);
  
      % And we can then define samplepoints for the fluorescence
      if (imoff-imon+1) >= stdall
          imstep = stdall;
      else
          imstep = imoff - imon + 1;
      end
      samplepoints = zeros(imstep,1);
      dt = floor((imoff - imon + 1)/imstep); 
      samplepoints(1) = imon;
      samplepoints(end) = imoff;
      
      for k = 2:imstep-1
          samplepoints(k) = round((k-1)*dt) + imon; 
      end
      
      % And... HOP! We extract corresponding fluorescence values and we
      % then interpolate data to the corresponding timestamps of the stride
      % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      isofluo(ss,:) = interp1(im.time(samplepoints),mintensity(samplepoints),timestamps(ss,:)); 
      isofluo(ss,:) = (isofluo(ss,:) - isofluo(ss,1))/isofluo(ss,1); 
      % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
  end
    
  isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated = isofluo; 
  isolated_fluorescence.([paw_dictionary{pw},'_paw']).timestamps = timestamps;  
end

%% V  Isolating fluorescence in groups of 'stdall' points for all cells  
% Now that we have extracted strides for each trial in the session, we'd
% like to snip the corresponding fluorescence values, either from bulk
% fluorescence or for individual cells. Identically, we will sample
% (stdall) points of fluorescence using the stride timestamps and store
% everything in a structure, 'isolated_fluorescence' with paw-specific
% subfields and 'isolated' + 'timestamps' further subsubfields. 

clear cell_isolated_fluorescence
Ncells = cnf.n_cells;

% just like for the previous struct 

% We will loop through paws. 
for pw = 1:4
  
  % First we extract the corresponding tracking data 
  isolated = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
  timestamps = isolated_strides.([paw_dictionary{pw},'_paw']).timestamps;
  
  %How many strides will we snip fluorescence for 
  [Nstr,~] = size(isolated);

  %preallocation
  isofluo = zeros(Nstr,Ncells,stdall); 
  spk = zeros(Nstr,Ncells,stdall); 
  imtimestamps = zeros(Nstr,stdall); 
  imbins = zeros(Nstr,1); 
  
  % We now loop through each stride of this paw 
  for ss = 1:Nstr
      onset = timestamps(ss,1); %refs(1)
      offset = timestamps(ss,end); %refs(3) 
      
      % In this paragraph, we find the closest timestep in im.time for both
      % the stride onset and offset.
      whereon = abs(im.time - onset);
      imon = find(whereon == min(whereon)) - 1; 
      imon = imon(1); 
      whereoff = abs(im.time - offset);
      imoff = find(whereoff == min(whereoff)) + 1; 
      imoff = imoff(1);
  
      % And we can then define samplepoints for the fluorescence
      if (imoff-imon+1) >= stdall
          imstep = stdall;
      else
          imstep = imoff - imon + 1;
      end
      samplepoints = zeros(imstep,1);
      dt = floor((imoff - imon + 1)/imstep); 
      samplepoints(1) = imon;
      samplepoints(end) = imoff;
      L = length(im.time(imon:imoff));
      
      
      for k = 2:imstep-1
          samplepoints(k) = round((k-1)*dt) + imon; 
      end
      
      for cc = 1:Ncells
          ivect = cnf.intensity(:,cc); 
          svect = cnf.spikes(:,cc);
          % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
          isofluo(ss,cc,:) = interp1(im.time(samplepoints),ivect(samplepoints),timestamps(ss,:));
          % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
          svect = svect(imon:imoff);
          propt = stdall/L;
          spkidx = find(svect == 1);
          svect = zeros(stdall,1);
          for ii = 1:length(spkidx)
              newidx = round(spkidx(ii) * propt);
              if newidx == 0
                  newidx = 1;
              end
              svect(newidx) = svect(newidx) + 1;
          end
          spk(ss,cc,:) = svect; 
          % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      end
      imbins(ss) = L;
  end
    
  cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated = isofluo; 
  cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).timestamps = timestamps; 
  cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).spikes = spk;
  cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).imbins = imbins;
end

%% VI  PLOTS 
%% ________________________________________________________________________
%  From there on we make use of the extracted data to analyse the
%  stride-locked trajectory of bulk fluorescence, cell-specific
%  fluorescence, spike probabilities, etc. 
%% a  Stride-locked bulk fluorescence
% What is the trajectory of mean bulk fluorescence values across the step
% cycle? 

figure, hold on
reordered = [3,1,4,2];

mvar = zeros(4,1); 
stvar = zeros(4,1);

for k = 1:4
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    subplot(2,2,k), hold on;
    fluodata = isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    [n,~] = size(fluodata);
    incertitude = std(fluodata,0,1)/sqrt(n); 
    strdata = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
    [Nstr,~] = size(fluodata);
    
     yyaxis right
     plot(mean(strdata,1),'LineWidth',2.5,'Color',ccc), hold on
     ax=gca; 
     ax.YAxis(1).Color = ccc;
     ylabel('x coord.','FontWeight','bold')
     
     
     yyaxis left
     %plot(mean(fluodata,1),'LineWidth',2.5,'Color',ccc)
     %errorbar(mean(fluodata,1),incertitude,'Color',ccc)
     clear LineProps
     LineProps.width = 1;
     LineProps.edgestyle = '-';
     LineProps.col = {ccc};
     mseb(linspace(1,10,10),mean(fluodata,1),incertitude,LineProps,1); hold on
     ax=gca; 
     ax.YAxis(2).Color = ccc;
     ax.XAxis.Color = ccc;
     ax.XTick = linspace(1,10,10);
     ax.XTickLabels = num2cell(linspace(1,10,10));
     ylabel('\Delta F/F') 
     mvar(pw,1) = mean(incertitude);
     stvar(pw,1) = std(incertitude)/sqrt(length(incertitude)); 
     ylim([-10e-3 5e-3])
    xlim([0 11])
    xlabel('Time bins')
    title([paw_dictionary{pw},' paw'],'Color',ccc)
    
    box off 
    
end

suptitle('Stride-locked bulk fluorescence') 

figure, hold on
box off
b = bar(linspace(1,4,4),mvar);
b.FaceColor = 'flat';
for k = 1:4
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    b.CData(pw,:) = ccc;
end
b.FaceAlpha = 0.5;
xticks(linspace(1,4,4))
xticklabels(paw_dictionary)
set(gca,'TickLength',[0,0.5])
xlabel('Paws')
ylabel('RMSE')
errorbar(mvar,stvar,'Color','k','LineStyle','none')
title('Amplitude of bulk-fluorescence modulation for the 4 paws')


%% b  Cell-resoluted mean fluorescence traces locked onto the stride cycle for each paw
% Here we prepare 4 matrices (1 for each paw) with individual cells in
% lines and the previously-used 10 time bins in columns containing the mean
% fluorescence value of the cell in question over all stride cycles. 

reordered = [3,1,4,2];
for k = 1:4
    
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    
    isolated = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    misolated = zeros(Ncells,stdall); 
    for cc = 1:Ncells
        cellspec = squeeze(isolated(:,cc,:)); 
        [nstr,~] = size(cellspec);
        ivect = zeros(nstr,stdall);
        maxint = max(cnf.intensity(:,cc));
        minint = min(cnf.intensity(:,cc));
        cellspec = (cellspec - minint)/ maxint;
        misolated(cc,:) = zero_and_max(mean(cellspec,1)); 
    end

    sub=subplot(2,2,k);
    c = gray(1000);
    c = c.*ccc; 
    colormap(gca,c)
    imagesc(misolated); ylabel('Cells'), title([paw_dictionary{pw},' paw'],'Color',ccc), %colorbar, box off,
    ax=gca; 
    ax.YAxis.Color = ccc;
    ax.XAxis.Color = ccc;
    hold off 
    
end

%% c  Reorganise cells based on metakmeans clusters
% We replot the 4 matrices but with clustered lines instead of monotonic
% ones. 

reordered = [3,1,4,2];
orderedcells = [];
K = length(activity_clustersK4.clusterregions);
oo = [3 2 1];
for i = 1:K
    k = oo(i);
    alldem = activity_clustersK4.clusterregions{1,k};
    orderedcells = cat(1,orderedcells,alldem); 
end

Nclcells = length(orderedcells);
for k = 1:4
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    
    isolated = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    misolated = zeros(Nclcells,stdall); 
    for cc = 1:Nclcells
        roi = orderedcells(cc);
        cellspec = squeeze(isolated(:,roi,:)); 
%         misolated(cc,:) = zero_and_max(mean(cellspec,1)); 
        [nstr,~] = size(cellspec);
        ivect = zeros(nstr,stdall);
        maxint = max(cnf.intensity(:,cc));
        minint = min(cnf.intensity(:,cc));
        cellspec = (cellspec - minint)/ maxint;
        misolated(cc,:) = zero_and_max(mean(cellspec,1)); 
    end

    sub=subplot(2,2,k);
    c = gray(1000);
    c = c.*ccc; 
    colormap(gca,c)
    imagesc(misolated); ylabel('Cells'), title([paw_dictionary{pw},' paw'],'Color',ccc), %colorbar, box off,
    ax=gca; 
    ax.YAxis.Color = ccc;
    ax.XAxis.Color = ccc;
    hold off 
    
end


clear orderedcells 
orderedcells = [];
for k = 1:3
    alldem = find(indx == k);
    orderedcells = cat(1,orderedcells,alldem); 
end

Nclcells = length(orderedcells);
for k = 1:4
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    
    isolated = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    misolated = zeros(K,stdall); 
    for i = 1:3
        alldem = find(indx == i);
        nn = length(alldem);
        subisolated = zeros(nn,stdall);
        for ii = 1:nn
            roi = alldem(ii);
            cellspec = squeeze(isolated(:,roi,:)); 
            subisolated(ii,:) =  zero_and_max(mean(cellspec,1));
        end
        misolated(i,:) = mean(subisolated,1);
    end

    sub=subplot(2,2,k);
    c = gray(1000);
    c = c.*ccc; 
    colormap(gca,c)
    imagesc(misolated); ylabel('Cluster'), title([paw_dictionary{pw},' paw'],'Color',ccc), %colorbar, box off,
    ax=gca; 
    ax.YAxis.Color = ccc;
    ax.XAxis.Color = ccc;
    hold off 
    
end

%% d    Separate cell fluorescence based on the metakmeans clusters
% We do so in shades of grey based on the previous meta-k-means clustering.
% We then compute a cluster-specific mean fluorescence trace to see if
% their locking to the stride differs. We finally re-plot the purkinje
% artscape with corresponding grey colours. 

figure, hold on
paw_colors = [0.976, 0.223, 0.243 ; 0.972, 0.266, 0.831 ; 0.223, 0.388, 0.976 ; 0.223, 0.866, 0.976];
reordered = [3,1,4,2];

time = isolated_fluorescence.([paw_dictionary{1},'_paw']).timestamps;
time = time(1,:); 
time = time - time(1); 
time = round(time,2);
clcolors = [0 0 0 ; 0.35 0.35 0.35; 0.6 0.6 0.6];

clmvar = zeros(4,3); 
clstvar = zeros(4,3); 

for k = 1:4
    
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    
    subplot(2,2,k)
    xlim([0,11])
    yyaxis right
    strdata = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
    plot(mean(strdata,1),'LineWidth',4,'Color',cat(2,ccc,0.45)), hold on
    ax=gca; 
    ax.YAxis(1).Color = ccc;
    ylabel('x coord.')
    
    isolated = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    misolated = zeros(Ncells,stdall); 
    for cc = 1:Ncells
        cellspec = squeeze(isolated(:,cc,:)); 
        mcellspec = mean(cellspec,1);
        misolated(cc,:) = (mcellspec - mcellspec(1))/mcellspec(1); 
    end

    
    yyaxis left
    for cc = 1:Ncells
        for j = 1:3
            cl = activity_clustersK4.clusterregions{1,j};
            %cl = find(indx == j);
            if ~isempty(find(cl==cc))
                clcolor = clcolors(j,:);
                clcolor = cat(2,clcolor,0.15);
                %plot(misolated(cc,:),'Color',clcolor,'Marker','none'), hold on 
            end
        end
    end
    
    for j = 1:3
            cl = activity_clustersK4.clusterregions{1,j};
            %cl = find(indx == j);
            clcolor = clcolors(j,:);
            %plot(mean(misolated(cl,:),1),'Color',clcolor,'Marker','none','LineWidth',1.5,'LineStyle','-'), hold on 
            incertitude = std(misolated(cl,:),0,1)/sqrt(length(cl)); 
            %errorbar(mean(misolated(cl,:),1),incertitude,'Color',clcolor,'Marker','none')
            clmvar(pw,j) = mean(incertitude);
            clstvar(pw,j) = std(incertitude)/sqrt(length(incertitude));
            clear LineProps
            LineProps.width = 1;
            LineProps.edgestyle = '-';
            LineProps.col = {clcolor};
            mseb(linspace(1,10,10),mean(misolated(cl,:),1),incertitude,LineProps,1); hold on
    end
    


    %plot(mean(misolated,1),'LineWidth',2.5,'Color',ccc)
    %errorbar(mean(misolated,1),incertitude,'Color',ccc)
    ax=gca; 
    ax.YAxis(2).Color = ccc;
    ax.XAxis.Color = ccc;
    ax.XTick = linspace(1,10,10);
    ax.XTickLabels = linspace(1,10,10);
    ylabel('\Delta F/F') 
    xlabel('Time bin')
    box off
    ylim([-12e-3 5e-3])

    hold off 
    
    title([paw_dictionary{pw},' paw'],'Color',ccc)
    
end

suptitle('Mean fluorescence traces for each of the 3 kmeans clusters')


figure, hold on
bposi = [1 2 3 ; 5 6 7 ; 9 10 11 ; 13 14 15];

for k = 1:4
    pw = reordered(k);
    ccc = paw_colors(pw,:);
    b = bar(bposi(pw,:),clmvar(pw,:));
    b.FaceColor = 'flat';
    for cl = 1:3
        b.CData(cl,:) = ccc + (cl-1) * 0.3;
    end
    b.FaceAlpha = 0.5;
    errorbar(bposi(pw,:),clmvar(pw,:),clstvar(pw,:),'Color','k','LineStyle','none');
end
xticks([2 6 10 14])
xticklabels(paw_dictionary)
xlabel('Paws')
ylabel('RMSE')


%% e   Random separation of cells in K groups (test of cluster significance)
% HERE WE RANDOMLY SPLIT CELLS IN 3 GROUPS TO SHOW THAT THE CLUSTERING DOES
% ACTUALLY GENERATE ASSEMBLIES OF PCS WITH A COMPACT BEHAVIOUR

figure, hold on
paw_colors = [0.976, 0.223, 0.243 ; 0.972, 0.266, 0.831 ; 0.223, 0.388, 0.976 ; 0.223, 0.866, 0.976];
reordered = [3,1,4,2];

time = isolated_fluorescence.([paw_dictionary{1},'_paw']).timestamps;
time = time(1,:); 
time = time - time(1); 
time = round(time,2);
clcolors = [0 0 0 ; 0.35 0.35 0.35; 0.6 0.6 0.6];

clmvar = zeros(4,3); 
clstvar = zeros(4,3); 


articlsize = [56 57 57];
population = linspace(1,170,170);
artifcl = cell(3,1);
for j = 1:3
    x = randsample(population,articlsize(j));
    artifcl(j,:) = {x};
    for k = 1:articlsize(j)
        whisit = find(population == x(k));
        population(whisit) = []; 
    end
end


for k = 1:4
    
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    
    subplot(2,2,k)
    xlim([0,11])
    yyaxis right
    strdata = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
    plot(mean(strdata,1),'LineWidth',4,'Color',cat(2,ccc,0.45)), hold on
    ax=gca; 
    ax.YAxis(1).Color = ccc;
    ylabel('x coord.')
    
    isolated = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    misolated = zeros(Ncells,stdall); 
    for cc = 1:Ncells
        cellspec = squeeze(isolated(:,cc,:)); 
        mcellspec = mean(cellspec,1);
        misolated(cc,:) = (mcellspec - mcellspec(1))/mcellspec(1); 
    end

    
    yyaxis left
  
    
    for j = 1:3
            cl = artifcl{j,:}; 
            %cl = find(indx == j);
            clcolor = clcolors(j,:);
            %plot(mean(misolated(cl,:),1),'Color',clcolor,'Marker','none','LineWidth',1.5,'LineStyle','-'), hold on 
            incertitude = std(misolated(cl,:),0,1)/sqrt(length(cl)); 
            %errorbar(mean(misolated(cl,:),1),incertitude,'Color',clcolor,'Marker','none')
            clmvar(pw,j) = mean(incertitude);
            clstvar(pw,j) = std(incertitude)/sqrt(length(incertitude));
            clear LineProps
            LineProps.width = 1;
            LineProps.edgestyle = '-';
            LineProps.col = {clcolor};
            mseb(linspace(1,10,10),mean(misolated(cl,:),1),incertitude,LineProps,1); hold on
    end
    


    %plot(mean(misolated,1),'LineWidth',2.5,'Color',ccc)
    %errorbar(mean(misolated,1),incertitude,'Color',ccc)
    ax=gca; 
    ax.YAxis(2).Color = ccc;
    ax.XAxis.Color = ccc;
    ax.XTick = linspace(1,10,10);
    ax.XTickLabels = linspace(1,10,10);
    ylabel('\Delta F/F') 
    xlabel('Time bin')
    box off
    ylim([-12e-3 5e-3])

    hold off 
    
    title([paw_dictionary{pw},' paw'],'Color',ccc)
    
end

suptitle('Mean fluorescence traces for each of the 3 kmeans clusters')


figure, hold on
bposi = [1 2 3 ; 5 6 7 ; 9 10 11 ; 13 14 15];

for k = 1:4
    pw = reordered(k);
    ccc = paw_colors(pw,:);
    b = bar(bposi(pw,:),clmvar(pw,:));
    b.FaceColor = 'flat';
    for cl = 1:3
        b.CData(cl,:) = ccc + (cl-1) * 0.3;
    end
    b.FaceAlpha = 0.5;
    errorbar(bposi(pw,:),clmvar(pw,:),clstvar(pw,:),'Color','k','LineStyle','none');
end
xticks([2 6 10 14])
xticklabels(paw_dictionary)
xlabel('Paws')
ylabel('RMSE')




%%
figure; hold on

title('Spatial distribution of clustered regions')
S = size(cnf.mask{1,1});
for k = 1:3
   regions = artifcl{k,:};
   L = length(regions);
   c = clcolors(k,:);
   full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
   for roi = 1:L
      actual = regions(roi);
      I = cnf.mask{1,actual};
      h = imshow(full); hold on 
      set(h, 'AlphaData', I*0.9) , hold on
   end
end

%% f  Does the variance of fluorescence traces follow any specific trend either? 
% In the initial bulk fluorescence x stride plot, we added error bars which
% amplitude did not seem to change much along the binning. We want to make
% sure of this and see how variance changes through the stride, meaning :
% do cells come closer to a stereotyped behaviour (synchrony?) at some
% times of the step cycle? This is one step closer to plotting spikes on
% top of this. 


reordered = [3,1,4,2];
for k = 1:4
    
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    sub=subplot(2,2,k);
    yyaxis right
    strdata = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
    plot(mean(strdata,1),'LineWidth',2,'Color',ccc), hold on
    ax=gca; 
    ax.YAxis(1).Color = ccc;
    ylabel('x coord.')
    
    yyaxis left
    isolated = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    [Nstr,~,~] = size(isolated); 
    stdisolated = zeros(Ncells,stdall); 
    for cc = 1:Ncells
        cellspec = squeeze(isolated(:,cc,:)); 
        for b = 1:stdall
            thebin = cellspec(:,b); 
            stdisolated(cc,b) = std(thebin)/sqrt(Nstr); %the variance for this cell and this bin
        end
    end
    
    mstdisolated = mean(stdisolated,1); %we average variances for bins across cells 
    incertitude = std(stdisolated,0,1)/sqrt(Ncells); %this is variance of variance! Hilarious. 
    c = gray(1000);
    c = c.*ccc; 
    colormap(gca,c)
    plot(mstdisolated,'Color',cat(2,ccc,0.5),'LineWidth',4); ylabel('\sigma'), title([paw_dictionary{pw},' paw'],'Color',ccc), %colorbar, box off,
    errorbar(mstdisolated,incertitude,'Color',ccc)
    ax=gca; 
    ax.YAxis(2).Color = ccc;
    ax.XAxis.Color = ccc;
    box off
    hold off 
    
    
end
suptitle('Evolution of fluorescence variance across all cells throughout the stride cycle')

%% g    Plot the mean fluorescence variation over all strides for a single cell 
% Here we just want to see to what extent a given cell has a stereotyped
% behaviour in its contribution to the observed stride modulation, meaning
% : is the trend that we see at the mean level already visible when we plot
% single-stride fluorescence traces, usw. The cell reference is specified
% in 'roi'. 

roi =  1;


figure, hold on
paw_colors = [0.976, 0.223, 0.243 ; 0.972, 0.266, 0.831 ; 0.223, 0.388, 0.976 ; 0.223, 0.866, 0.976];
reordered = [3,1,4,2];

time = isolated_fluorescence.([paw_dictionary{1},'_paw']).timestamps;
time = time(1,:); 
time = time - time(1); 
time = round(time,2);

for k = 1:4
    
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    
    subplot(2,2,k)
    xlim([0,11])
    yyaxis right
    strdata = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
    plot(mean(strdata,1),'LineWidth',4,'Color',cat(2,ccc,0.45)), hold on
    ax=gca; 
    ax.YAxis(1).Color = ccc;
    ylabel('x coord.')
    
    isolated = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    cellspec = squeeze(isolated(:,roi,:));
    [nstr,stdall] = size(cellspec); 
    misolated = mean(cellspec,1); 
   

    
    yyaxis left
    
    for s = 1:nstr %use this to plot all stride-locked traces 
        %plot(cellspec(s,:),'Color',[0.6 0.6 0.6 0.5]), hold on
    end
    
    plot(misolated,'LineWidth',2.5,'Color',ccc)
    errorbar(misolated,std(cellspec,0,1)/sqrt(nstr),'Color',ccc)
    ax=gca; 
    ax.YAxis(2).Color = ccc;
    ax.XAxis.Color = ccc;
    ax.XTick = linspace(1,10,10);
    ax.XTickLabels = num2cell(time);
    ylabel('\Delta F/F') 
    xlabel('Time (s)')
    box off

    hold off 
    

    title([paw_dictionary{pw},' paw'],'Color',ccc)
    
end

%% h     Amplitude of the stride-locked fluorescence modulation
% Here, I calculate the mean stride modulation for each cell for each paw's
% cycles and I define the fluorescence modulation, which is the peak (up)
% to peak (down) y distance (normalised fluorescence). We then see if
% modulation is greater in a paw, in a side of the body and if all cells
% seem to contribute equally (gaussian vs bimodal vs ...). We then fit the
% distribution and do a bar plot of the mean modulation amplitude. t-tests
% and ANOVA are then used to test the statistical significance of the
% observed differences. 


reordered = [3,1,4,2];
clcolors = [0 0 0 ; 0.35 0.35 0.35; 0.6 0.6 0.6];

pkdistrib = zeros(4,Ncells); %we use this to store the modulation values
%for each cell and each paw 

% - - - - - - - - - - - - - - HISTOGRAMS - - - - - - - - - - - - - - - - - 
for k = 1:4
    
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    
    pktopk = zeros(Ncells,1);
    isolated = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    for cc = 1:Ncells
        % . . . . . . With this each cell has a proper normalisation . . . 
        cellspec = squeeze(isolated(:,cc,:));
        cellspec = (cellspec - min(cnf.intensity(:,cc)))/max(cnf.intensity(:,cc));
        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        misolated = mean(cellspec,1); 
        pktopk(cc) = max(misolated) - min(misolated);  
    end

    probdist = fitdist(pktopk,'Normal'); %we fit the dist. with a normal law (after visual inspection) 
    xvals = 0:0.00001:0.015;
    yfit = pdf(probdist,xvals);
    
    subplot(2,2,k)
    
    % PLOTTING
    bw = 0.00025;
    yyaxis left
    h = histogram(pktopk,'Normalization','probability','Binwidth',bw); hold on
    [N,~] = histcounts(pktopk,'Normalization','pdf','Binwidth',bw);
    whereismax = find(N == max(N));
    mode = whereismax(1) * bw;
    h.FaceColor = ccc;
    h.FaceAlpha = 0.5; 
    h.EdgeColor = ccc - 0.2; 
    ax=gca; 
    ax.YAxis(1).Color = ccc;
    ax.XAxis.Color = ccc;
    xlim([0,0.015])
    yyaxis right
    plot(xvals,yfit,'Color',cat(2,ccc,0.6),'LineWidth',1.8)
    ax.YAxis(2).Color = ccc;
    ax.YTick(:) = [];
    ylabel('Probability') 
    xlabel(sprintf('Normalised amplitude \n (of the mean stride modulation)'))
    box off

    hold off 
    
    pkdistrib(k,:) = pktopk; %and store the paw-specific distribution of amps. 
    title([paw_dictionary{pw},' paw - mode = ',num2str(mode),],'Color',ccc)
    
end

% - - - - - - - - - - - - - - BARPLOT - - - - - - - - - - - - - - - - - - - 
figure, hold on
for k = 1:4
    pw = reordered(k);
    ccc = paw_colors(pw,:);
    b = bar(k,mean(pkdistrib(k,:)));
    b.FaceColor = ccc;
    b.FaceAlpha = 0.6;
    b.EdgeColor = ccc; - 0.1;
    errorbar(k,mean(pkdistrib(k,:)),std(pkdistrib(k,:))/sqrt(Ncells),'Color',ccc)
end
xticks(1:4)
xticklabels(paw_dictionary(reordered))
ylabel('Mean stride-locked amplitude modulation (\Delta F/F)')
xlabel('Paws')

%% i   Spike histograms (PSTH-like)
% Here we use the deconvolution data to see where spike localise during the
% step cycle.

figure, hold on
paw_colors = [0.976, 0.223, 0.243 ; 0.972, 0.266, 0.831 ; 0.223, 0.388, 0.976 ; 0.223, 0.866, 0.976];
reordered = [3,1,4,2];

time = isolated_fluorescence.([paw_dictionary{1},'_paw']).timestamps;
time = time(1,:); 
time = time - time(1); 
time = round(time,2);
allprobspikes = zeros(4,Ncells,stdall);

for k = 1:4
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    subplot(2,2,k);
    spikes = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).spikes;
    strdata = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
    imbins = cell_isolated_fluorescence.([paw_dictionary{pw},'_paw']).imbins;
    Nimbins = sum(imbins);
    [Nstr,Ncells,stdall] = size(spikes);
    
     yyaxis right
     plot(mean(strdata,1),'LineWidth',2,'Color',cat(2,ccc,0.6)), hold on
     ax=gca; 
     ax.YAxis(1).Color = ccc;
     ylabel('x coord.','FontWeight','bold')
     
     
     yyaxis left
     probspikes = zeros(Ncells,stdall);
     for ii = 1:Ncells
         allspikes = cnf.spikes(:,ii);
         Pspk = length(find(allspikes == 1))/length(allspikes); %baseline spiking probability
         
         for bb = 1:stdall
             binspk = squeeze(spikes(:,ii,bb)); 
             Pbinspk = sum(binspk)*10/Nimbins;
             probspikes(ii,bb) = Pbinspk/Pspk;
             %probspikes(ii,bb) = Pbinspk;
         end
     end
     
     mprobspikes = mean(probspikes,1);
     b = bar(mprobspikes);
     incertitude = std(probspikes,0,1)/sqrt(Ncells); 
     errorbar(mprobspikes,incertitude,'Color',ccc,'LineStyle','none')
     b.FaceColor = ccc;
     b.FaceAlpha = 0.7;
     b.EdgeColor = ccc - 0.1; 
     ax=gca; 
     ax.YAxis(2).Color = ccc;
     ax.XAxis.Color = ccc;
     ax.XTick = linspace(1,10,10);
     ax.XTickLabels = num2cell(linspace(1,10,10));
     ylabel('P_{spike}^{norm}') 

     
    xlim([0 11])
    xlabel('Time bins')
    title([paw_dictionary{pw},' paw'],'Color',ccc)
    
    box off 
    allprobspikes(k,:,:) = probspikes;
end


for k = 1:4
    figure, hold on
    pw = reordered(k);
    ccc = paw_colors(pw,:); 
    probspikes = squeeze(allprobspikes(k,:,:));
    
    for bin = 1:stdall
        subplot(2,5,bin)
        h = histogram(probspikes(:,bin),'Normalization','count','Binwidth',0.05);
        title(['Bin ',num2str(bin)],'Color',ccc-0.1)
        h.FaceColor = ccc; 
        h.FaceAlpha = 0.6;
        h.EdgeColor = ccc - 0.1;
        ax = gca;
        ax.XAxis.Color = ccc;
        ax.YAxis.Color = ccc;
        box off 
        if bin == 1 || bin == 6
            ylabel('Number of cells')
            xlabel('Probability modulation')
        end
    end
    suptitle([paw_dictionary{pw},' paw'])
    
   
    hold off 

end
%%
pw_strides = [161,85,176,102];
b = bar(pw_strides);
b.FaceColor = 'flat';
for k = 1:4
    b.CData(k,:) = paw_colors(k,:);
end
b.FaceAlpha = 0.7;
box off 
xlabel('Paws')
xticklabels(paw_dictionary)
ylabel('Number of strides')