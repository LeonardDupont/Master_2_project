%% Giving a go to the MLSpike algorithm (Deneux et al. 2016)


%%
fs = 30; %Hz

%if we approximate that we can use the same kernel for all regions of
%interest (not likely), then we can loop the autocalibrations on ROIs
taus = zeros(20,1);
amps = zeros(20,1);

for roi=1:20
    data = calcium_data(:,roi);
    [tau, amp, sigma, events] = spk_autocalibration(data,1/fs);
    taus(roi) = tau;
    amps(roi) = amp; 
    clear tau, clear amp
end

tau_m = mean(taus); 
amp_m = mean(amps); 

%%
% and then we deconvolute the calcium traces
par = tps_mlspikes('par')
par.dt = 1/30;
par.a = 
par.tau = 
par.saturation =
par.finetune.sigma = 
par.drift.parameter = 


for roi=1:20 
    [spk, fit, drift, parest] = spk_est(calcium_data(:,roi),par); 
    deconvolution.(['roi_',num2str(roi)]).spiketimes = spk;  %we store everything in a big struct with rois as subfields
    deconvolution.(['roi_',num2str(roi)]).fit = fit;
    %deconvolution.(['roi_',num2str(roi)]).drift = drift; 
    deconvolution.(['roi_',num2str(roi)]).sigma = parest.finetune.sigma; 
end

%% 
% visualisation
roi = 1; %choose the region (can eventually be looped) 
spike_t = deconvolution.(['roi_',num2str(roi)]).spiketimes; 
spk_display(par.dt,spike_t,calcium_data(:,roi))