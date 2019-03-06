function [cn] = denoise_and_demix(cn,neuropile,cn_bkg,varargin)
% March 2019 - Carey lab - leonard.dupont@ens.fr
% .........................................................................
% This function is the final step in our attempt at demixing and denoising
% the fluorescence signals output by the Mukamel ROIs during 1P-imaging.
% At this stage, in the cn struct, the fluorescence intensity Fi of cell i 
% at frame (time) t can be expressed as explained hereunder: 
%
%
%           Fi(t) = [ki*si](t) + ci*Ni(t) + sigma(t)
%
%           where:      ki     is the calcium transient kernel
%                       si     is the spiketrain (convolution w\ ki)
%                       Ni     is the neuropile fluorescence (mixed)
%                       sigma  is the noise
%
%   <=>     [ki*si](t) = Fi(t) - ci*Ni(t) - sigma(t)
%
% Here we basically demix and denoise the roi's signal at each frame by 
% implementing the above-written substraction. Indeed, we have a background
% roi and we just wrote functions to extract a neuropile fluorescence for 
% each putative purkinje unit.
% .........................................................................
%
% -- INPUT -------
%
%     cn            struct with fields (contains Fi)
%     neuropile     struct with fields (contains Ni)
%     cn_bkg        struct with fields (contains sigma)
%     varargin      ci value, set to 0.7 by default
%
% -- OUTPUT -------
% 
%     cn           struct with identical fields, but added [ki*si](t)
%                  (cn.intensity_dmdn : 'demixeddenoised')
% .........................................................................

ip = inputParser;
ip.addParameter('ci',0.7); %maybe we'll find a way to optimise this
parse(ip,varargin{:})
ci = ip.Results.ci;


%% Building cn.intensity_dmdn
N = cn.n_cells;
frames = length(cn.intensity);

sigma = cn_bkg.intensity;
for roi = 1:N
    Fi = cn.intensity{1,roi};
    Ni = neuropile.intensity{1,roi};
    intensity_dmdn = zeros(frames,1);
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    for t = 1:frames
        intensity_dmdn(t,1) = Fi(t) - ci*Ni(t) - sigma(t);
    end
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  
    cn.intensity_dmdn{1,roi} = intensity_dmdn;

end