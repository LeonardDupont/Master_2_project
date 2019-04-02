function [speedkk] = imbased_konkat(im,tm,N)
%% April 2019 - CareyLab - leonard.dupont@ens.fr
% .........................................................................
% This function uses the discontinuities in the im timestamps to detect
% trial begins, trial ends and concatenated the speed data in groups of 
% N trials (N should be the same as in concatenate_mukamel.m). 
% .........................................................................
%
%    INPUT
%
%    im       struct with fields as created by Hugo's functions (imaging)
%    tm       same but for treadmill
%     N       number of trials in each contatenation group (MATCH MUKAMEL!)
%
%    OUTPUT
%
%   speedkk  speed(KonKatenated), either a vector (nG = 1) or a struct
%
% .........................................................................

if nargin < 3
    N = 4;
    warning('No concatenation size input by user. Using 4 by default...')
end

%% 1 - finds discontinuities and marks the times
startpts = [];
stoppts = [];

fs = 30; %Hz
epsilon = 2;
startpts(end+1) = im.time(1);

for k = 1:length(im.time)-1
    delta = im.time(k+1) - im.time(k);
    if (delta - 1/fs) > epsilon
        startpts(end+1) = im.time(k+1);
        stoppts(end+1) = im.time(k);
    end
end
stoppts(end+1) = im.time(end);

startpts= startpts* 1/tm.delta_time;
stoppts= stoppts* 1/tm.delta_time; %conversion to frames (speed) using fs2

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

if nG > 1 %if more than one group, then we will output a struct
    
    for k = 1:nG
        concat_speed = [];
        if k == 1
            initial = 1;
            final = groupsize(1);
        else
            initial = groupsize(k-1) + 1;
            final = final + groupsize(k);
        end
        for j = initial:final
            start = floor(startpts(j));
            stop = floor(stoppts(j)); 
            speedd = tm.speedM(start:stop);
            concat_speed = cat(1,concat_speed,speedd); 
        end
    end
    speedkk.(['kk',num2str(k)]) = concat_speed;
    
else %if  there is only one group
    speedkk = [];
    initial = 1;
    final = N;
    for j = initial:final
            start = floor(startpts(j));
            stop = floor(stoppts(j)); 
            speedd = tm.speedM(start:stop);
            speedkk = cat(1,speedkk,speedd); 
    end
end

disp('Done!')
end