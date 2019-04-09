function [ISIfreq, average_spkrate] = build_ISI_histo(data,varargin)
% This function takes in spiketimes and builds the interspike-interval
% histogram. February 2019 - Carey lab (LD)
%
% INPUT
% 
% data            can either be a vector (single trace spiketimes) or the
%                 cn struct with field 'spikes'
%
% varargin        see below
% 
%
% OUTPUT
%
% the function plots a histogram of ISI distribution probabilities


fs = 30; %30 Hz default 
ip = inputParser; 
ip.addParameter('fs',fs, @isscalar); %the framerate if different from the default one 
ip.addParameter('individual',0); %if you feed the function with the struct but you want to calculate separate histograms for each ROI
ip.addParameter('graphics',1); %if you want to plot the histogram or not (might no want to if individual is 1)
ip.addParameter('bw',0.1, @double); %the bin width in the histogram
ip.addParameter('color','k')
ip.addParameter('rmhigh',[])
parse(ip, varargin{:});

fs = ip.Results.fs;
individual = logical(ip.Results.individual);
graphics = logical(ip.Results.graphics);
bw = ip.Results.bw;
colour = ip.Results.color;
rmhigh = ip.Results.rmhigh;

%% 
if isstruct(data) 
    
    N = data.n_cells;
    cn = data;
    clear data
    
    if individual
        
        for roi=1:N
            
            spikes = cn.spikes(:,roi); 
            spike_times = find(spikes == 1);
            ISIfreq = zeros(length(spike_times)-1,1);
            
            for t=1:length(spike_times)-1
                ISIfreq(t) = fs/(spike_times(t+1) - spike_times(t));
            end
            
                            if ~isempty(rmhigh)
                                artefacts = find(ISIfreq > rmhigh);
                                ISIfreq(artefacts) = [];
                            end
            
            if graphics
                figure; hold on
                title(['ISI histogram for ROI ', num2str(roi)])
                histogram(ISIfreq,'BinWidth',bw)
                 xlabel('ms bin')
            end
            
            [N,edges] = histcounts(ISIfreq,'BinWidth',bw,'FaceColor',colour,'Normalization','probability');
            mm = median(ISIfreq); 
        end
        
        
    else %if we can use them all 
        ISIfreq = [];
        
        for roi=1:N
            spikes = cn.spikes(:,roi); 
            spike_times = find(spikes == 1);
            for t=1:length(spike_times)-1
                ISI = spike_times(t+1) - spike_times(t);
                ISIfreq = [ISIfreq, fs/ISI]; 
            end
            
                            if ~isempty(rmhigh)
                                artefacts = find(ISIfreq > rmhigh);
                                ISIfreq(artefacts) = [];
                            end 
            
        end
        if graphics
            figure; hold on
            title('ISI histogram all ROIs')
            histogram(ISIfreq,'BinWidth',bw,'FaceColor',colour,'Normalization','probability')
             xlabel('Hz bin')
        end
        
        [N,edges] = histcounts(ISIfreq,'BinWidth',bw);
        mm = median(ISIfreq); 
        
    end
       
else %then it's just a vector
    

    spike_times = find(data == 1); 
    ISIfreq = zeros(length(spike_times)-1,1);
    for t = 1:length(spike_times)-1
        ISIfreq(t) = fs/(spike_times(t+1)-spike_times(t));
    end
    
                            if ~isempty(rmhigh)
                                artefacts = find(ISIfreq > rmhigh);
                                ISIfreq(artefacts) = [];
                            end
    
    
    if graphics
            figure; hold on
            title('ISI histogram')
            histogram(ISIfreq,'BinWidth',bw,'FaceColor',colour,'Normalization','probability')
            xlabel('Hz bin')
    end   
    
    [N,edges] = histcounts(ISIfreq,'BinWidth',bw);
    mm = median(ISIfreq); 
    
end

max_bin = find(N == max(N));
average_spkrate = edges(max_bin); %we find the class
disp(['Max bin is for ',num2str(average_spkrate), ' Hz'])
disp(['Median frequency is ',num2str(mm),'Hz.'])

end


