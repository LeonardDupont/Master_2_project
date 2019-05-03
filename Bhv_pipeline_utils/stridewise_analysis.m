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
