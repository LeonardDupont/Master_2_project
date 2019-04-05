

% b : now we feed the autocalibration and deconvolute the signals
clear par
[x,y] = size(calcium_data);
good_trace = cn.intensity(:,57);

% autocalibration
calcium_norm = good_trace/mean(good_trace); %doesn't necessarily need to be normalised

pax = spk_autocalibration('par');
pax.dt = 1/30; %framerate
pax.realspikes = sspike_events*1/30; %the transients we deem as representative of sspike events (see above)
pax.mlspikepar.dographsummary = false;

[tauest, aest, sigmaest] = spk_autocalibration(calcium_norm,pax);

% deconvolution using Viterbi algorithm
par = tps_mlspikes('par');
par.a = aest;
par.dt = 1/30;
par.tau = tauest;
par.drift.parameter = 0.001;
p2 = 0.5;
p3 = 0.01;
par.pnonlin = [p2 p3] ;
par.dographsummary = false;
par.finetune.sigma = sigmaest; %to be left empty from my experience, even
%though we estimated it 

cn = rmfield(cn,'spikes');
for roi=1:cn.n_cells
    disp(['Processing ROI ' , num2str(roi)])
    [spk, fit, drift, parest] = spk_est(cn.intensity(:,roi),par); 
    cn.spikes{1,roi} = spk;
end