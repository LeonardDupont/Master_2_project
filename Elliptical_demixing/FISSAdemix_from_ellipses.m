function [cn] = FISSAdemix_from_ellipses(neuropile,cn)
%% March 2019 - CareyLab - leonard.dupont@ens.fr
% .........................................................................
% This function takes in the neuropile struct containing cell-specific info
% about the surrounding contamination and uses NNMF with NNDSVD
% initialisation to factorise the mixing signals from the center ROI.
% Inspired from Keeminsk et al. (2018).
% .........................................................................
%
%   INPUT
%
%   neuropile      struct with fields
%   cn             struct with fields
%
%
%   OUTPUT
%
%   cn             identical struct with added field 'intensity_dm'. 
%
% .........................................................................

nseg = neuropile.nseg;
if neuropile.bkg
    nseg = nseg+1;
end


N = neuropile.n_cells;
% . . . . . . . Normalise before demixing, otherwise artefacts . . . . . . 
for roi = 1:N
    for seg = 1:(nseg+1)
        a = zero_and_max(neuropile.intensity{seg,roi}); 
        neuropile.intensity{seg,roi} = a ; 
    end
end



% . . . . . . . . NNMF occurs here with NNDSVD initialisation . . . . . . .
for roi = 1:N
    
   if rem(roi,10) == 0 || roi == N
       disp(['Cell ',num2str(roi),' out of ',num2str(N),'.'])
   end
    
    tic
    disp('1 -- Preparing fluorescence matrix with neuropile subregions --')
    F = initFISSA_mixedF(neuropile,roi);
    toc
    disp('2 -- NNDSVD-based initialisation --')
    [W0,H0] = NNDSVD(F,nseg+1,0); 
    %initialisation of weight and separated matrices using SVD 
    

    % . . . . . . . . minimising |F - W*H| according to Frobenius . . . . .
    [W,H] = nnmf(F,nseg+1,'w0',W0,'h0',H0); 
   
    if neuropile.bkg
        maxi = max(W(nseg,:)); %then last but one segment (last = bkg)
        wc = find(W(nseg,:) == maxi);
    else
        maxi = max(W(nseg+1,:)); %then last segment
        wc = find(W(nseg+1,:) == maxi);
    end
    dmxd = H(wc,:);
    
    cn.intensity_dm(:,roi) = dmxd; 
    
end




end