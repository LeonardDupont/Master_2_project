%% APPLYING CNMF-E POST MUKAMEL - February 2019, Carey lab -- leonard.dupont@ens.fr
%
%   Part 1 : creating CNMF-E - compatible objects 
%   Part 2 : CNMF-E by Zhou et al. (2018)
%   Part 3 : transforming the objects back to Carey neurons
%
%% 1 - building initialisation inputs from the mukamel algorithm
% Right after running Mukamel, we have a 'cn' and a 'cn_bkg' struct. 
% The CNMF-E algorithm needs the following inputs to run:
%
%   A       d * K matrix (spatial components, unwrapped frames)
%   b       d * 1 matrix (background component, unwrapped frame)
%   C       K * T matrix (temporal components)
%   f       T * 1 matrix (background temporal info)
%   center  K * 2 matrix (centroid coordinates {x,y} forall K)


filepath = '/Users/leonarddupont/Desktop/M2_internship/[Clustering] Figures_and_docs/MC318_26_218_F_tied_138,000_137,000_5_1_m_Ready.tif';
Y = imread_tifflib(filepath); %fast
[d1,d2,d3] = size(Y);
%Y = Y - min(Y(:));
T = d3;

%% a. preprocess data : get parameters
p = 2; %order of the autoregressive model
[Y,P] = preprocess_data(Y,p); %seems necessary to get P, can also take a look at its format and build it myself

%% b. building a sparse matrix (Ain) for spatial components
K = cn.n_cells;
d = d1*d2;
A = zeros(d,K); 

for roi = 1:K
    S = cn.roi_landscape{1,roi};
    S = reshape(S,[d 1]); %unwrap it to make it a (d1*d2)-long vector 
    A(:,roi) = S;
end

Ain = sparse(A); %for storage purposes 

%% c. building (Cin) and (center) for temporal components and ROI {x,y}

Cin = zeros(K,T); %temporal components : fluorescent trace
center = zeros(K,2); %centroids

for roi = 1:K
    Cin(roi,:) = (cn.intensity(:,roi)).';
    centroid = cn.centroid{1,roi};
    center(roi,1) = centroid(1,1);
    center(roi,2) = centroid(1,2);
end

%% d. estimating background components (in space and time) using NMF

Yr = reshape(Y,numel(Y)/T,T); %now a d * T matrix
medY = median(Y,2); %median of the tiff stack over time 
Y = bsxfun(@minus, Y, medY); %we remove the median 
nb = 1; %1 bkg component


[bin,fin] = nmf(max(Y - Ain*Cin + repmat(medY,1,T),0),nb); 



%% 2 - Updating spatial and temporal components iteratively to demix signals
% We have now transformed the Mukamel-extracted data into CNMF-E
% compatible objects. We will now iteratively apply the algorithm on our
% ROIs to demix signals from spatial neighbours. For each repeat:
%
%   (i)     We update the spatial components
%   (ii)    We update the temporal components
%   (iii)   Overlapping components are ROIs and are merged
%

repeats = 3;

A = Ain;
C = Cin;
b = bin;
f = fin;


for run = 1:repeats
    
    P.p = p;
    [A,b,C] = update_spatial_components(Yr,C,f,A,P,options);
    
    if run == 1
        P.p = 0;    % apparently just for speed : quicky first run
        [C,f,P,S] = update_temporal_components(Yr,A,b,C,f,P,options);
    else
        [C,f,P,S] = update_temporal_components(Yr,A,b,C,f,P,options);
    end
    
    [Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);
    A = Am;
    C = Cm;
    K = K_m; %preparing the next iteration  
    
end



%% 3 - We put the data back in the cn format
% Now we can transform the CNMF-E - processed data back into Carey neurons.

cn_dmxd.fov_width = w;
cn_dmxd.fov_height = h;
cn_dmxd.n_cells = K; %straightforward and common to all rois


Af = full(A); 
Af = reshape(Af,[d1,d2,d3]); %back to the frame shape

for roi = 1:K
    cn_dmxd.roi{1,:} = find(A(:,roi) ~= 0); %just locations of non-0 pixel
    cn_dmxd.roi_landscape{1,roi} = Af(:,:,roi); %a full frame with values
    
    premask = Af(:,:,roi);
    premask(premask ~= 0) = 1; %replacing the landscape values with 1s
    cn_dmxd.mask{1,:} = premask;
    
    cn_dmxd.intensity = C(roi,:); %raw trace
    
    [x,y] = ind2sub([d1 d2],cn_dmxd.roi{1,:});
    xc = sum(x)/length(x);
    yc = sum(y)/length(y); %calculating the centroid's x and y coordinates in the frame
    cn_dmxd.centroid{1,roi}(1,1) = xc;
    cn_dmxd.centroid{1,roi}(1,2) = yc;
end





