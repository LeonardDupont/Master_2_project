%% Processing session data

session_raw_dir = 'Z:\leonard.dupont\TM RAW FILES\voluntary locomotion\MC318\S5';
tracking_dir = 'Z:\leonard.dupont\TM TRACKED FILES\MC318\S5';
imaging_dir = 'Z:\leonard.dupont\TM IMAGING FILES\voluntary locomotion\MC318\S5';
session_output_dir = 'Z:\leonard.dupont\TM SESSION FILES\voluntary locomotion\';

[exp_files] = get_experimental_files_ordered_by_animal_session_and_trial(session_raw_dir, '*.tdms');
platform = get_default_widefield_rotary_treadmill_parameters(2);

concatenate_session_TDMS_data(exp_files, session_output_dir, platform);
concatenate_session_RTM_treadmill_data(session_output_dir, platform);
concatenate_session_RTM_tracking_data(exp_files, tracking_dir, session_output_dir, platform);
concatenate_session_RTM_imaging_data(exp_files, imaging_dir, session_output_dir, platform);

load('/Volumes/carey/leonard.dupont/TM SESSION FILES/voluntary locomotion/MC318/S5/imaging_data.mat')
load('/Volumes/carey/leonard.dupont/TM SESSION FILES/voluntary locomotion/MC318/S5/port_data.mat')
load('/Volumes/carey/leonard.dupont/TM SESSION FILES/voluntary locomotion/MC318/S5/tracking_data.mat')
load('/Volumes/carey/leonard.dupont/TM SESSION FILES/voluntary locomotion/MC318/S5/treadmill_data.mat')

%% Initial RTM analysis (26/03/2019)
calcium_data = cn.intensity.';
m_ca = mean(calcium_data,1);
m_ca = zero_and_max(m_ca);

% Initial plot : over the whole session 
figure, hold on
plot(tm.time,tm.speedM), axis tight
plot(im.time,m_ca)

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
speedkk = imbased_konkat(im,tm,4); 
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

tm_speedMs = smoothdata(speedkk,'gaussian',10000); %smoothed version
figure, plot(tm_speedMs)

spdtime = linspace(1,tm.time(length(tm_speedMs)),length(tm_speedMs)); 
imtime = linspace(1,tm.time(length(tm_speedMs)),length(m_ca));

rescale = max(tm_speedMs);

figure, hold on
plot(spdtime,tm_speedMs), axis tight
plot(imtime,m_ca*rescale)


% . . . . . . . . SPEED CATEGORIES ANALYSIS . . . . . . . . . . . . . . . . 
ncat = 3; 

speedcats = define_speedcats(tm_speedMs,ncat);

speedlabels = assign_speedlabels(tm_speedMs,speedcats);
ncat = max(speedlabels); %taking in2 account the 0s (make up 1 cat)

speedlabels = get_speedregime_boundaries(speedlabels);


figure, hold on
clear cell 
clear l 
cmap = jet(ncat);
legend_spd = cell(ncat,1);
for i = 1:ncat
    l(i) = plot([NaN,NaN], 'color', cmap(i,:));
    spdlbl = ['speed_',num2str(i-1)];
    legend_spd(i) = cellstr(spdlbl); 
end
legend(l, legend_spd,'AutoUpdate','off');

for k = 1:length(xc)-1
    xplus = xc(k+1);
    xminus = xc(k);
    label = speedlabels(xplus-1);
    color = cmap(label,:);
    
    xl = [xminus, xminus, xplus, xplus];
    yl = [0 100 100 0];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.4;
    a.EdgeColor = 'none';
end
        
axis tight, hold on
plot(linspace(1,length(speedlabels),length(tm_speedMs)),tm_speedMs*(100/max(tm_speedMs)),':','color','k','LineWidth',1)


plot(linspace(1,length(speedlabels),length(m_ca)),m_ca*(100/max(m_ca)),'color','k','LineWidth',0.001)
box off 
yticklabels = [];

%% Mutual information analysis
Tim = length(m_ca);
Ttm = length(tm_speedMs);
mergeF = Ttm/Tim; 

binboundaries = zeros(k,1); 
for k = 1:Tim
    binboundaries(k) = floor((k-1)*mergeF + 1);
end

aside = Ttm - binboundaries(end);
binboundaries(end+1) = binboundaries(end)+aside;
figure, hold on
Nclust = length(activity_cluster.clusterregions);
for j = 1:Nclust
    regions = activity_cluster.clusterregions{:,j};
    N = length(regions);
    allclusttraces = zeros(length(cn.intensity(:,1)),N);
    for roi = 1:N
        allclusttraces(:,roi) = zero_and_max(cn.intensity(:,regions(roi)).').';
    end
    mclusttrace = mean(allclusttraces,2);
    mclusttrace = mclusttrace.';
    %for roi = 1:N
        %trace = cn.intensity(:,roi).';
        %trace = zero_and_max(trace);
        fluocats = define_fluocats(mclusttrace,10);
        fluolabels = assign_fluolabels(mclusttrace,fluocats);

        Nimbins = max(fluolabels);
        Ntmbins = max(speedlabels);
        MImat = zeros(Nimbins,Ntmbins);

        for k = 1:Tim
            start = binboundaries(k);
            stop = binboundaries(k+1);
            tmlabel = round(mean(speedlabels(start:stop)));
            imlabel = fluolabels(k);
            MImat(imlabel,tmlabel) = MImat(imlabel,tmlabel) + 1;
        end

        MImat = MImat/Tim; 

        for imlbl = 1:Nimbins
            Pim = length(find(fluolabels == imlbl))/ Tim;
            for tmlbl = 1:Ntmbins
                Ptm = length(find(speedlabels == tmlbl))/ Ttm;
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                MImat(imlbl,tmlbl) =  log(MImat(imlbl,tmlbl)/ (Pim*Ptm));
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            end
        end





    subplot(2,3,j),imagesc(MImat)
    title(['Cluster ',num2str(j)])
    ylabel('Fluorescence bins')
    xlabel('Speedbins')
    if rem(j,Nclust/2) == 0
        colorbar
    end
end
suptitle('Point mutual information in all clusters')

%%


