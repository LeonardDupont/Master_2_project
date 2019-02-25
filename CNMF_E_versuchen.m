%% 1 - building initialisation inputs from the mukamel algorithm
% Right after running mukamel, we have a 'cn' and a 'cn_bkg' struct. 
% The CNMF-E algorithm needs the following inputs to run:
%
%   A       d * K matrix (spatial components)
%   b       d * 1 matrix (background component)
%   C       K * T matrix (temporal components)
%   f       T * 1 matrix (background temporal info)
%   center  K * 2 matrix (centroid coordinates {x,y} forall K)
% 
% We will build this from our structs. 

filepath = '/Users/leonarddupont/Desktop/M2_internship/[Clustering] Figures_and_docs/MC318_26_218_F_tied_138,000_137,000_5_1_m_Ready.tif';
Y = imread_tifflib(filepath);
Y = Y - min(Y(:));
T = ndim(Y);

%% preprocess data : get parameters
p = 2; %order of the autoregressive model
[Y,P] = preprocess_data(Y,p);

%% we must now build a sparse matrix Ain corresponsing to the ROI landscapes
% we them extracted using Mukamel (2009). 
K = cn.n_cells;
d1 = cn.fov_height;
d2 = cn.fov_width;
d = d1*d2;

A = zeros(d,K);
for roi = 1:K
    S = cn.roi_landscape{1,roi};
    S = reshape(S,[d 1]); 
    A(:,roi) = S;
end

Ain = sparse(A);
%% we can now also build Cin and center, which is far simpler

Cin = zeros(K,T);
center = zeros(K,2);
for roi = 1:K
    Cin(roi,:) = (cn.intensity(:,roi)).';
    centroid = cn.centroid{1,roi};
    center(roi,1) = centroid(1,1);
    center(roi,2) = centroid(1,2);
end

%% Final initialisation : estimation of background 
Yr = reshape(Y,numel(Y)/T,T); %now a d * T matrix
medY = median(Y,2); %median of the tiff stack over time 
Y = bsxfun(@minus, Y, medY); %we remove the median 
nb = 1; %1 background matrix


[bin,fin] = nmf(max(Y - Ain*Cin + repmat(medY,1,T),0),nb); 


%% 2 - Updating spatial and temporal components iteratively to demix signals
% We have now transformed the Mukamel-extracted data into a CNMF-E
% compatible one.
repeats = 3;

A = Ain;
C = Cin;
b = bin;
f = fin;


for run = 1:repeats
    P.p = p;
    [A,b,C] = update_spatial_components(Yr,C,f,A,P,options);
    if run == 1
        P.p = 0;    % set AR temporarily to zero for speed
        [C,f,P,S] = update_temporal_components(Yr,A,b,C,f,P,options);
    else
        [C,f,P,S] = update_temporal_components(Yr,A,b,C,f,P,options);
    end
    
    [Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);
    
    A = Am;
    C = Cm;
    K = K_m;
   
end

