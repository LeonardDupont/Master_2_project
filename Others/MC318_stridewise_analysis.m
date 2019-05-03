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

tracking_dir = 'Z:\leonard.dupont\TM PROCESSED FILES\voluntary locomotion'; 
tdmsdata = 'Z:\leonard.dupont\TM SESSION FILES\voluntary locomotion\MC318\S4'; 
tm = load(fullfile(tdmsdata,'tracking_data.mat')); 
tm = tm.tracks; 
im = load(fullfile(tdmsdata,'imaging_data.mat'));
im = im.im; 
filePattern = fullfile(tracking_dir,'*.mat');
Ntrials = length(dir(filePattern)); 
trials = dir(filePattern); 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

%%
fstm = 300; %Hz
fsim = 30; %Hz

tmidx = find_trial_onset(tm.time,fstm); 
imidx = find_trial_onset(im.time,fsim);

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
stdpts1 = 4; 
stdpts2 = 3;
stdall = stdpts1 + stdpts2 + 3; 

countit = [0 0 0 0]; 

clear isolated_strides
paw_dictionary = {'FR','HR','FL','HL'}; 


for k = 1:Ntrials
    realidx = find(order == k); 
    name = trials(realidx).name; 
    disp(['[stride extraction] Processing ', name, '.'])
    dtMAT = load([tracking_dir,'\',name]);
    
    strides = dtMAT.strides; 
    paws = dtMAT.trial_data.paws_x;
    
    if k < Ntrials
        ttime = tm.time(tmidx(k):tmidx(k+1)); 
    else
        ttime = tm.time(tmidx(k):end);
    end
    
    for pw = 1:4
        coord = strides{pw}; 
        Nstr = length(coord); 
        isolated = zeros(Nstr,stdall);
        timestamps = zeros(Nstr,stdall); 
        
        for st = 1:Nstr

            refs = coord(st,:); 
            
            samplepoints = zeros(stdall,1);
            samplepoints(1) = refs(1); 
            samplepoints(end) = refs(3);
            samplepoints(1 + stdpts1 + 1) = refs(2); %refs(1) + stdpts + 1
            
            dt1 = (refs(2)-refs(1))/(stdpts1+1); 
            dt2 = (refs(3)-refs(2))/(stdpts2+1); 
            for dd = 1:stdpts1
                samplepoints(1+dd) = round(dd*dt1) + refs(1); 
            end 
            
            for dd = 1:stdpts2
                samplepoints(dd+stdpts1+2) = round(dd*dt2) + refs(2);
            end

            isolated(st,:) = paws(samplepoints,pw);
            timestamps(st,:) = ttime(samplepoints);
        end
        
        if countit(pw)
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

%%
clear isolated_fluorescence
countit = [0 0 0 0]; 

for pw = 1:4
  isolated = isolated_strides.([paw_dictionary{pw},'_paw']).isolated;
  timestamps = isolated_strides.([paw_dictionary{pw},'_paw']).timestamps;
  [Nstr,~] = size(data);

  isofluo = zeros(Nstr,stdall); 
  
  for ss = 1:Nstr
      onset = timestamps(ss,1); 
      offset = timestamps(ss,end); 
      
      whereon = abs(im.time - onset);
      imon = find(whereon == min(whereon)) - 1; 
      whereoff = abs(im.time - offset);
      imoff = find(whereoff == min(whereoff)) + 1; 
  
      dt = (imoff - imon)/stdall; 
      samplepoints = zeros(stdall,1);
      for k = 1:stdall
          samplepoints(k) = round((k-1)*dt) + imon; 
      end
      
      isofluo(ss,:) = intensity(samplepoints); 
      isofluo(ss,:) = interp1(isofluo(ss,:),im.time(samplepoints),timestamps); 
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



