%autocalibration
calcium_norm = good_trace/mean(good_trace);

pax = spk_autocalibration('par');
pax.dt = 1/30;
pax.amin = 0.003;
pax.amax = 1;
pax.taumin = 0.1;
pax.taumax = 1.8;
pax.saturation = 7e-4; %GCaMP6f
pax.realspikes = sspike_events;
pax.reala = [];
pax.mlspikepar.dographsummary = false;

[tauest, aest, sigmaest] = spk_autocalibration(calcium_norm,pax);
%%

pcal = spk_calibration('par');
pcal.dt = 1/30 ;
pcal.tdrift = 1; 
pcal.dosaturation = 7e-4;
pcal.dodelay = 15.6e-3;
pcal.dohill = 2.99;
pcal.doc0 = 1.5;
spikes = selected;

[pest, fit, drift, rcalcium] = spk_calibration(spikes,calcium_norm,pcal);


%%
clear par
par = tps_mlspikes('par');
par.a = 0.15;
par.dt = 1/30;
par.tau = tauest;
par.drift.parameter = 0.0001;
par.dographsummary = false;

[spikest, fit, drift, parest] = spk_est(good_trace,par); 

figure(1)
spk_display(1/30,{spikest},{good_trace fit drift})
set(1,'numbertitle','off','name','Autocal + MLspike')

%%

for roi=1:20 
    [spk, fit, drift, parest] = spk_est(calcium_data(:,roi),par); 
    deconvolution.(['roi_',num2str(roi)]).spiketimes = spk;  %we store everything in a big struct with rois as subfields
    deconvolution.(['roi_',num2str(roi)]).fit = fit;
    %deconvolution.(['roi_',num2str(roi)]).drift = drift; 
    deconvolution.(['roi_',num2str(roi)]).sigma = parest.finetune.sigma; 
end

%% visualisation
roi = 3; %choose the region (can eventually be looped) 
spike_t = deconvolution.(['roi_',num2str(roi)]).spiketimes; 
spk_display(par.dt,spike_t,calcium_data(:,roi))

%%
[pks,locs] = define_single_events(very_good);
%%
sspike_events = select_single_spikes(locs);

