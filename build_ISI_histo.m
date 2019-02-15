function [interspike, average_spkrate] = build_ISI_histo(data,varargin)
% This function takes in spiketimes and builds the interspike-interval
% histogram. February 2019 - Carey lab (LD)
%
% INPUT
% 
% data            can either be a vector (single trace spiketimes) or a matrix
%                 in which case the histogram will be more statistically
%                 relevant (the deconvolution matrix)
%
% varargin        see below
% 
%
% OUTPUT
%
% the function plots a histogram of ISI distribution probabilities


dt = 1/30; %30 Hz default 
ip = inputParser; 
ip.addParameter('dt',dt, @isscalar); %the framerate if different from the default one 
ip.addParameter('individual',0); %if you feed the function with the struct but you want to calculate separate histograms for each ROI
ip.addParameter('graphics',1); %if you want to plot the histogram or not (might no want to if individual is 1)
ip.addParameter('bw',0.1, @double); %the bin width in the histogram
ip.addParameter('color','k')
ip.addParameter('rmhigh',[])
parse(ip, varargin{:});

dt = ip.Results.dt;
individual = logical(ip.Results.individual);
graphics = logical(ip.Results.graphics);
bw = ip.Results.bw;
colour = ip.Results.color;
rmhigh = ip.Results.rmhigh;

%% 
if isstruct(data) %then it is a deconvolution-like struct 
    
    n_rois = numel(fieldnames(data));
    
    if individual
        
        for roi=1:n_rois
            
            spike_times = data.(['roi_',num2str(roi)]).spiketimes;
            interspike = zeros(length(spike_times)-1,1);
            
            for t=1:length(spike_times)-1
                interspike(t) = 1/(spike_times(t+1) - spike_times(t));
            end
            
                            if ~isempty(rmhigh)
                                artefacts = find(interspike > rmhigh);
                                interspike(artefacts) = [];
                            end
            
            if graphics
                figure; hold on
                title(['ISI histogram for ROI ', num2str(roi)])
                histogram(interspike,'BinWidth',bw)
                 xlabel('Hz bin')
            end
            
            [N,edges] = histcounts(interspike,'BinWidth',bw,'FaceColor',colour,'Normalization','probability');
            
        end
        
        
    else %if we can use them all 
        interspike = [];
        
        for roi=1:n_rois
            spike_times = data.(['roi_',num2str(roi)]).spiketimes;
            for t=1:length(spike_times)-1
                ISI = 1/(spike_times(t+1) - spike_times(t));
                interspike = [interspike, ISI]; 
            end
            
                            if ~isempty(rmhigh)
                                artefacts = find(interspike > rmhigh);
                                interspike(artefacts) = [];
                            end 
            
        end
        if graphics
            figure; hold on
            title('ISI histogram all ROIs')
            histogram(interspike,'BinWidth',bw,'FaceColor',colour,'Normalization','probability')
             xlabel('Hz bin')
        end
        
        [N,edges] = histcounts(interspike,'BinWidth',bw);
        
    end
       
else
    
    interspike = zeros(length(data)-1,1);
    for t = 1:length(data)-1
        interspike(t) = 1/(data(t+1)-data(t));
    end
    
                            if ~isempty(rmhigh)
                                artefacts = find(interspike > rmhigh);
                                interspike(artefacts) = [];
                            end
    
    
    if graphics
            figure; hold on
            title('ISI histogram')
            histogram(interspike,'BinWidth',bw,'FaceColor',colour,'Normalization','probability')
            xlabel('Hz bin')
    end   
    
    [N,edges] = histcounts(interspike,'BinWidth',bw);
    
end

max_bin = find(N == max(N));
average_spkrate = edges(max_bin); %we find the class
disp(['Max bin is for ',num2str(average_spkrate), ' Hz'])

end


