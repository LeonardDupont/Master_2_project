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
% First, we build a cn-like struct with fluorescence data from the whole
% session. Fluorescence is extracted using the usual masks. 

wholeS4 = bigread2(''); %import tiff 
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


%% Extract strides

% loading imaging, tracking tdms as well as tracking files 
tracking_dir = 'Z:\leonard.dupont\TM PROCESSED FILES\voluntary locomotion'; 
tdmsdata = 'Z:\leonard.dupont\TM SESSION FILES\voluntary locomotion\MC318\S4'; 
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

fstm = 300; %Hz
fsim = 30; %Hz

tmidx = find_trial_onset(tracks.time,fstm); 
imidx = find_trial_onset(im.time,fsim);

%% Isolating strides in groups of 'stdall' points
% Here, we travel through all the session trials and extract stride points
% (x values) as well as corresponding timestamps in the session. Everything
% is stored in a struct 'isolated_strides' with subfields for different
% paws. For each paw, we then have a matrix 'isolated' with datapoints of
% every detected stride for the limb in the session and a matrix 
% 'timestamps' with corresponding time values. 


% To sample the step cycle, we use three reference points : stance offset (1),
% swing (2) and stance onset (3). We add 4 points between (1) and (2) and 3
% points between (2) and (3). This can be changed hereunder of course. 
stdpts1 = 4; 
stdpts2 = 3;
stdall = stdpts1 + stdpts2 + 3; 

countit = [0 0 0 0]; %this is to correctly build the struct (see below)
% we count if this the first time we treat the paw in question 

clear isolated_strides
% 1 : Front Right -- 2 : Hind Right -- 3 : Front Left -- 4 : Hind Left
paw_dictionary = {'FR','HR','FL','HL'}; 

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
            samplepoints(1 + stdpts1 + 1) = refs(2); %refs(1) + stdpts + 1

            %first bit
            dt1 = (refs(2)-refs(1))/(stdpts1+1); %refs(1) to refs(2)
            for dd = 1:stdpts1
                samplepoints(1+dd) = round(dd*dt1) + refs(1); 
            end 
            
            %second bit
            dt2 = (refs(3)-refs(2))/(stdpts2+1); %refs(2) to refs(3)
            for dd = 1:stdpts2
                samplepoints(dd+stdpts1+2) = round(dd*dt2) + refs(2);
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

%%
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

%% Isolating fluorescence in groups of 'stdall' points 
% Now that we have extracted strides for each trial in the session, we'd
% like to snip the corresponding fluorescence values, either from bulk
% fluorescence or for individual cells. Identically, we will sample
% (stdall) points of fluorescence using the stride timestamps and store
% everything in a structure, 'isolated_fluorescence' with paw-specific
% subfields and 'isolated' + 'timestamps' further subsubfields. 

clear isolated_fluorescence

% just like for the previous struct 
countit = [0 0 0 0]; 

% We will loop through paws. 
for pw = 1:4
  
  % First we extract the corresponding tracking data 
  isolated = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
  timestamps = isolated_strides.([paw_dictionary{pw},'_paw']).timestamps;
  
  %How many strides will we snip fluorescence for 
  [Nstr,~] = size(data);

  %preallocation
  isofluo = zeros(Nstr,stdall); 
  
  % We now loop through each stride of this paw 
  for ss = 1:Nstr
      onset = timestamps(ss,1); %refs(1)
      offset = timestamps(ss,end); %refs(3) 
      
      % In this paragraph, we find the closest timestep in im.time for both
      % the stride onset and offset.
      whereon = abs(im.time - onset);
      imon = find(whereon == min(whereon)) - 1; 
      whereoff = abs(im.time - offset);
      imoff = find(whereoff == min(whereoff)) + 1; 
  
      % And we can then define samplepoints for the fluorescence 
      dt = (imoff - imon)/stdall; 
      samplepoints = zeros(stdall,1);
      for k = 1:stdall
          samplepoints(k) = round((k-1)*dt) + imon; 
      end
      
      % And... HOP! We extract corresponding fluorescence values and we
      % then interpolate data to the corresponding timestamps of the stride
      % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      isofluo(ss,:) = intensity(samplepoints); 
      isofluo(ss,:) = interp1(isofluo(ss,:),im.time(samplepoints),timestamps); 
      % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
  end
    
  imtimestamps = im.time(samplepoints); 
  if countit(pw)
     isofluo = cat(1,isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated,isofluo); 
     imtimestamps = cat(1,isolated_fluorescence.([paw_dictionary{pw},'_paw']).timestamps,imtimestamps); 
  end
  isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated = isofluo; 
  isolated_fluorescence.([paw_dictionary{pw},'_paw']).timestamps = imtimestamps;  
  countit(pw) = 1; 
end


%% Plots 

figure, hold on

for pw = 1:4
    subplot(2,2,pw)
    fluodata = isolated_fluorescence.([paw_dictionary{pw},'_paw']).isolated;
    strdata = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
    [Nstr,~] = size(fluodata);
    
    for k = 1:Nstr
        yyaxis right
        plot(strdata(k,:),'Color',[0.8 0.8 0.8]), hold on
        yyaxis left 
        plot(fluodata(k,:),'Color',[0.6 0 0]), hold on
    end
    
    yyaxis right
    plot(mean(strdata,1),'LineWidth',1.5,'Color','k'), hold on
    ylabel('x coord.')
    yyaxis left 
    plot(mean(fluodata,1),'LineWidth',1.5,'Color',[0.9 0 0])
    ylabel('\Delta F / F') 
    
    title([paw_dictionary{pw},' paw'])
end

suptitle('Stride-locked bulk fluorescence') 

%% Then we plot that for all cells, do a PSTH, and so on. Nice stuff 






