function [grdfov] = grid_grouping(data,gridbins)
%% March 2019 - Carey lab - leonard.dupont@ens.fr
%..........................................................................
% This function takes the movie in data and applies a grid with dimensions 
% requested by the user to divide the FOV in squares of fused pixels. 
% ........................................................................
% 
% -- INPUT ------------
%
%   data        h * w * time, tiff movie
%
%   gridbins   double, number of squares in x and y directions (identical)
%
%
% -- OUTPUT ------------    
%
%   grdfov   struct with fields gridpx.intensity, gridpx.centroid
%
% ........................................................................
disp(['Grouping pixels in ',num2str(gridbins*gridbins),' squares in the FOV.'])
tic
%% First defining boundary coordinates of each grid square
[h,w,t] = size(data);

total = gridbins*gridbins;

hmerged = floor(h/gridbins);
haside = rem(h,gridbins);
hcoord = zeros(gridbins,1);

wmerged = floor(w/gridbins);
waside = rem(w,gridbins);
wcoord = zeros(gridbins,1);

for k = 1:gridbins
    hcoord(k) = k*hmerged + (k==gridbins)*haside;
    wcoord(k) = k*wmerged + (k==gridbins)*waside;
end

%% Centroids
centroid = cell(1,total);
c = 1;
for i = 1:gridbins
    for j = 1:gridbins
        if i == 1
            height = hcoord(i)/2;
        else
            height = (hcoord(i-1) + hcoord(i))/2;
        end
        if j == 1
            width = wcoord(j)/2;
        else
            width = (wcoord(j-1) + wcoord(j))/2;
        end
        centroid{1,c} = [width , height]; %reads in height and then w
        c = c+1;
    end
end

%% Fluorescence intensity

intensity = zeros(total,t);

for frame = 1:t
    c = 1;
    for i = 1:gridbins
        for j = 1:gridbins
            
            if i == 1
               hstart = 1;
               hstop = hcoord(i);
            else
               hstart = hcoord(i-1) + 1;
               hstop = hcoord(i);
            end
            
            if j == 1
                wstart = 1;
                wstop = wcoord(j);
            else
                wstart = wcoord(j-1) + 1;
                wstop = wcoord(j);
            end
            
            
            thebox = data( hstart:hstop , wstart:wstop ,frame);
            intensity(c,frame) = mean(thebox(:)); 
            c = c+1;
        end
    end  
end

grdfov.centroid = centroid;
grdfov.intensity = intensity;

disp('Done!') 
toc

end