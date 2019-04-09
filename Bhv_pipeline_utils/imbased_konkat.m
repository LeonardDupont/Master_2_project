function speedkk = imbased_konkat(im,tm,inputpath,N)
%% April 2019 - CareyLab - leonard.dupont@ens.fr
% .........................................................................
% This function uses the discontinuities in the im timestamps to detect
% trial begins, trial ends and concatenated the speed data in groups of 
% N trials (N should be the same as in concatenate_mukamel.m). 
% CAREFUL, ONLY WORKS WITH TRIAL SIZE = 30 seconds
% Concatenation order is the same as in mukamel. 
% .........................................................................
%
%    INPUT
%
%    im       struct with fields as created by Hugo's functions (imaging)
%    tm       same but for treadmill
%     N       number of trials in each contatenation group (MATCH MUKAMEL!)
% inputpath   path to folder enclosing the '*Ready.tif' video files 
%
%    OUTPUT
%
%   speedkk  speed(KonKatenated), either a vector (nG = 1) or a struct
%
% .........................................................................

if nargin < 4
    N = 4;
    warning('No concatenation size input by user. Using 4 by default...')
end

trsize = 30; %in seconds, the duration of one trial in the behavioural paradigm


%% 1 - finds discontinuities and marks the times
startpts = [];
stoppts = [];

fs = 30; %Hz
epsilon = 2;
startpts(end+1) = im.time(1);
c=1;
for k = 1:length(im.time)-1
    delta = im.time(k+1) - im.time(k);
    if (delta - 1/fs) > epsilon
        startpts(end+1) = im.time(k+1);
        stoppts(end+1) = im.time(k);
        c = k;
    else
        if (k - c) == (trsize * 30)
            startpts(end+1) = im.time(k+1);
            stoppts(end+1) = im.time(k);
        end
    end
end
stoppts(end+1) = im.time(end);

startpts= startpts* 1/tm.delta_time;
stoppts= stoppts* 1/tm.delta_time; %conversion to frames (speed) using fs2

%% Building a list of trial numbers (for some are not saved and skipped)

filePattern = fullfile(inputpath, '*Ready.tif');

V = length(dir(filePattern));
if V == 0
    error('No registered tif files in this directory')
end
videofiles = dir(filePattern); 
trialnb = zeros(V,1);

for k = 1:V
    thename = videofiles(k).name;
    underscore = find(thename == '_');
    trial_nb = thename(underscore(end-2)+1:underscore(end-1)-1); 
    trialnb(k) = str2double(trial_nb);
end

%% 2 - Building groups of size N while accounting for remainer

L = length(startpts);
if N > L
    warning('Concatenation group size greater than number of trials.')
    disp(['Using N = ',num2str(L),' instead.'])
    N = L;
    groupsize = N;
else
    aside = rem(L,N);
    % this here is really to match what could have happened in
    % concatenate_mukamel.m (see function for details). 
    if aside < 4
        groupsize = N * linspace(1,1,floor(L/N));
        groupsize(end) = groupsize(end) + aside;
    else
        groupsize = N * linspace(1,1,floor(L/N));
        groupsize(end+1) = aside;
    end
end


%% 3 - Now concatenating speed-value trials to match the imaging data

nG = length(groupsize);

    clear speedkk
    for k = 1:nG
        concat_speed = [];
        if k == 1
            initial = 1;
            final = groupsize(1);
        else
            initial = 1 + sum(groupsize(1:k-1));
            final = sum(groupsize(1:k));
        end
        for j = initial:final
            start = floor(startpts(j));
            stop = floor(stoppts(j)); 
            speedd = tm.speedM(start:stop);
            concat_speed = cat(1,concat_speed,speedd); 
        end
    whichtrials = trialnb(initial:final);
    namest = 'trials';
    for i = 1:length(whichtrials)
        namest = [namest,'_',num2str(whichtrials(i))];
    end
    speedkk.(namest) = concat_speed;
    end
    
    speedkk.trialnb = trialnb;

disp('Done!')
end