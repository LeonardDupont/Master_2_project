
% again here, we must define single events for the autocalibration 
roi = 45; %check with the plot_all_traces function
good_trace = cn.intensity(:,roi);
[pks,locs] = define_single_events(good_trace);

% execute separately - GetGlobal command
 

% b : now we feed the autocalibration and deconvolute the signals 
clear deconvolution
clear par
[x,y] = size(calcium_data);


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


for roi=1:185
    disp(['Processing ROI ' , num2str(roi)])
    [spk, fit, drift, parest] = spk_est(calcium_data(:,roi),par); 
    cn.spiketimes{1,roi} = spk;
end