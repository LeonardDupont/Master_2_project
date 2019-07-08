load('port_data.mat')
port_datas4 = port_data ; 
clear port_data

load('treadmill_data.mat')
tm4 = tm; 
clear tm 


load('imaging_data.mat')
im4 = im;
clear im 


%%
figure, plot(tm4.time,tm_speedM), hold on
plot(im.time,cnf.intensity(:,1)*1000)
%%
fsim = 30; %Hz
imidx = find_trial_onset(im4.time,fsim);
imidx(end+1) = length(im4.time);
%%
tm_speedM = tm4.speedM;
ss_speedM = [];

for k = 1:length(imidx)-1
    imon = im4.time(imidx(k));
    imoff = im4.time(imidx(k+1)-1);
    
    closeon = abs(tm4.time - imon);
    closeoff = abs(tm4.time - imoff);
    whereon = find(closeon == min(closeon));
    whereoff = find(closeoff == min(closeoff));
    
    isolatedtrial = tm_speedM(whereon(1):whereoff(1));
    ss_speedM = cat(1,ss_speedM,isolatedtrial); 
    
end

%%
ss_speedM(ss_speedM < 0) = 0;
ss_speedMs = smoothdata(ss_speedM,'gaussian',2000);
%%
ncat = 1;
speedcats = define_speedcats(ss_speedMs,ncat);
speedlabels = assign_speedlabels(ss_speedMs,speedcats);
xc = get_speedregime_boundaries(speedlabels);


plot_all_traces_w_speed(cnf.intensity(:,1:30),xc,speedlabels,ss_speedMs)
%plot_spikes_w_speed(cnf.spikes,xc,speedlabels,ss_speedM)
%%
[~,Ncells] = size(cnf.intensity);
% Two different ways of computing activity in the 2 states, w/ different
% results and meanings. 
imlabelbounds = [];
imboundaries = [];
% - - - - - - - - - - - - - - - -
Allspkloco = [];
Nspkloco = cell(Ncells,1);
Allspkrest = [];
Nspkrest = cell(Ncells,1);
% - - - - - - - - - - - - - - - -
Pspkloco = zeros(Ncells,1); 
Pspkrest = zeros(Ncells,1); 
% - - - - - - - - - - - - - - - -

ttime = tm4.time(1:length(ss_speedM));
imtime = linspace(1,length(cnf.intensity(:,1))*1/30,length(cnf.intensity(:,1)));

Nlocobins = 0;
Nrestbins = 0;

% First we need to see where xc boundaries are in the imaging trace
for k = 1:length(xc)-1
    xminus = xc(k);
    xplus = xc(k+1);
    label = speedlabels(xplus-1);
    
    tmon = ttime(xminus);
    tmoff = ttime(xplus);
    closeon = abs(imtime-tmon);
    closeoff = abs(imtime-tmoff);
    
%     whereon = find(closeon == min(closeon));
%     whereoff = find(closeoff == min(closeoff));
    whereon = round(xminus/100)+1;
    whereoff = round(xplus/100);
    Nb = whereoff-whereon+1;
    if label == 1
        Nrestbins = Nrestbins + Nb;
        for roi = 1:Ncells
            spikes = cnf.spikes(whereon:whereoff,roi);
            whspk = find(spikes==1) + whereon - 1; 
            Nspk = length(whspk);
            Pspkrest(roi) = Pspkrest(roi) + Nspk;
            for i = 1:Nspk
                coactivated = sum(cnf.spikes(whspk(i),:)) - 1;
                Nspkrest{roi,:} = cat(2,Nspkrest{roi,:},coactivated); 
            end
        end
        for bb = 1:Nb
            Allspkrest(end+1) = sum(cnf.spikes(whereon + bb - 1,:));
        end
    else
        Nlocobins = Nlocobins + whereoff-whereon+1;
        for roi = 1:Ncells
            spikes = cnf.spikes(whereon:whereoff,roi);
            whspk = find(spikes==1) + whereon - 1; 
            Nspk = length(whspk);
            Pspkloco(roi) = Pspkloco(roi) + Nspk;
            for i = 1:Nspk
                coactivated = sum(cnf.spikes(whspk(i),:)) - 1;
                Nspkloco{roi,:} = cat(2,Nspkloco{roi,:},coactivated); 
            end
        end
        for bb = 1:Nb
            Allspkloco(end+1) = sum(cnf.spikes(whereon + bb - 1,:));
        end
    end
end

Pspkloco = Pspkloco / Nlocobins;
Pspkrest = Pspkrest / Nrestbins;


%%

% probability of spiking
figure, hold on
subplot(1,2,1)
h = histogram(Pspkrest,'Normalization','probability','Binwidth',0.003);
h.FaceColor = colors(1,:);
h.FaceAlpha = 0.6; 
h.EdgeColor = colors(1,:); hold on
h = histogram(Pspkloco,'Normalization','probability','Binwidth',0.003);
h.FaceColor = colors(2,:);
h.FaceAlpha = 0.6;
h.EdgeColor = colors(2,:);
title('PDist of spike likelihood')
box off 
xlim([0.01 0.095])
ylabel('Probability')
xlabel('P(spike)')

% subplot(1,3,2)
% h = histogram(Pspkloco,'Normalization','pdf','Binwidth',0.003);
% h.FaceColor = colors(2,:);
% h.FaceAlpha = 0.5;
% title('Moving','color',colors(2,:))
% xlabel('P(spike)')
% box off

subplot(1,2,2)
meanprobs = [mean(Pspkrest),mean(Pspkloco)];
incertitude = [std(Pspkrest)/sqrt(Ncells),std(Pspkloco)/sqrt(Ncells)];
b = bar(meanprobs);hold on
for roi = 1:Ncells
    plot([1,2],[Pspkrest(roi),Pspkloco(roi)],'Color',[0.8 0.8 0.8 0.5]), hold on
    scatter([1,2],[Pspkrest(roi),Pspkloco(roi)],5,[0.5 0.5 0.5],'filled')
end
errorbar(meanprobs,incertitude,'LineStyle','None','Color','k'), 
b.FaceColor = 'flat';
b.FaceAlpha = 0.6;
b.CData(1,:) = colors(1,:);
b.CData(2,:) = colors(2,:);
b.EdgeColor = 'none';
b.BarWidth = 0.9;
ylabel('P(spike)')
box off 
xticklabels({'Resting','Moving'})
set(gca,'TickLength',[0,0.01])
title('Overall spike probabilities')

suptitle('Impact of locomotion state on complex-spike probability')

%% Synchrony plot


figure, hold on
xposi = [1.5 4.5 7.5 10.5 13.5 16.5];
xbars = [1 2 ; 4 5 ; 7 8 ; 10 11 ; 13 14 ; 16 17];
xnames = {'1-10','10-20','20-30','30-40','50-60','>60'};
probchecks = [1 10 ; 10 20 ; 20 30 ;  30 40 ; 50 60 ; 60 170];
[l,~] = size(probchecks); 

for k = 1:l
    probdist = zeros(cnf.n_cells,2);
    nsynch = probchecks(k,:);
    for roi = 1:cnf.n_cells
            synchrest = Nspkrest{roi,:};
            interS = intersect(find(synchrest >= nsynch(1)),find(synchrest <= nsynch(2))); 
            probdist(roi,1) = length(interS)/length(synchrest);
            synchloco = Nspkloco{roi,:};
            interS = intersect(find(synchloco >= nsynch(1)),find(synchloco <= nsynch(2)));
            probdist(roi,2) = length(interS)/length(synchloco);
    end
  
    mrest = mean(probdist(:,1));
    mloco = mean(probdist(:,2)); 
    incertitude = [std(probdist(:,1))/sqrt(cnf.n_cells) , std(probdist(:,2))/sqrt(cnf.n_cells)];
    b = bar(xbars(k,:),[mrest,mloco]); 
    errorbar(xbars(k,:),[mrest,mloco],incertitude,'Color','k','LineStyle','none')
    b.FaceColor = 'flat';
    b.FaceAlpha = 0.6;
    b.BarWidth = 0.9;
    b.CData(1,:) = colors(1,:);
    b.CData(2,:) = colors(2,:);
    b.EdgeColor = 'none'; hold on
end

xticks(xposi)
xticklabels(xnames)
set(gca,'TickLength',[0, 0.05])
xlabel('Number of synchronous dendrites')
ylabel('Probability')
box off


figure
for k = 1:l
    probdist = zeros(cnf.n_cells,2);
    nsynch = probchecks(k,:);
    for roi = 1:cnf.n_cells
            synchrest = Nspkrest{roi,:};
            interS = intersect(find(synchrest >= nsynch(1)),find(synchrest <= nsynch(2))); 
            probdist(roi,1) = length(interS)/length(synchrest);
            synchloco = Nspkloco{roi,:};
            interS = intersect(find(synchloco >= nsynch(1)),find(synchloco <= nsynch(2)));
            probdist(roi,2) = length(interS)/length(synchloco);
    end
    subplot(2,3,k), hold on
    h1 = histogram(probdist(:,1),'Normalization','Probability');
    h1.FaceColor = colors(1,:);
    h1.FaceAlpha = 0.5;
    h1.NumBins = 15;
    h2 = histogram(probdist(:,2),'Normalization','Probability');
    h2.FaceColor = colors(2,:);
    h2.FaceAlpha = 0.5;
    h2.NumBins = 10;
    [p,h] = ranksum(probdist(:,1),probdist(:,2)); 
    title([xnames(k),'-',num2str(p)])
    box off 
end

%% With the first bissectrix
figure, hold on
x = linspace(0,0.1,100);
plot(x,x,'Color',[0.8 0.8 0.8])
for k = 1:cnf.n_cells
    scatter(Pspkrest(k),Pspkloco(k),15,[0 0 0],'filled'), hold on
end
axis tight
xlim([0 0.1])
ylim([0 0.1])
box off
xlabel('CS probability at rest')
ylabel('CS probability when moving')
xticks([0 0.02 0.04 0.06 0.08 0.1])
xticklabels({'0','','','','',0.1})
yticks([0 0.02 0.04 0.06 0.08 0.1])
yticklabels({'0','','','','',0.1})

%% Acceleration, decceleration and constant speed
% Here we define three acceleration regimes and see how these relate to
% cell activity. 

accs = get_acceleration_from_speed(ss_speedMs,tm4.time(1:length(ss_speedMs)));
mm = mean(accs);
accs = accs - mm;  %center
sigmaup = 0.02;
sigmadown = -0.02;

figure, plot(accs,'color','k','LineWidth',0.9), axis tight, hold on
yticks([])
xticks([])
box off 
xlabel('Time')
ylabel('Acceleration')
hold off 

%%
% Now we classify data 
acclabels = zeros(length(accs),1);
for k = 1:length(accs)
    pt = accs(k);
    if (pt<sigmaup) && (pt>sigmadown)
        acclabels(k) = 2; %small fluctuations
    elseif pt > sigmaup
        acclabels(k) = 3; %noteworthy acceleration
    elseif pt < sigmadown
        acclabels(k) = 1; %noteworthy decceleration
    end
end

acccolors = zeros(length(accs),3);
for k = 1:length(acclabels)
    switch acclabels(k)
        case 1
            acccolors(k,:) = [0 0.3 0.9];
        case 3
            acccolors(k,:) = [1 0.4 0.3];
        case 2
            acccolors(k,:) = [0.8 0.8 0.8];
    end
end
xcacc = get_speedregime_boundaries(acclabels);
fluo_x_accelerationpeaks(cnf.intensity,accs,acccolors,xcacc)

%% PMI mean trace
m_ca = mean(cnf.intensity,2); 

naccat = 3;
colormap = hot(naccat);

ttime = tm4.time(1:length(ss_speedM)-1);
imtime = linspace(0,ttime(end),length(m_ca));
rsm_ca = interp1(imtime,m_ca,ttime);

nfluocats = 5;
fluocats = define_fluocats(rsm_ca,nfluocats);
fluolabels = assign_fluolabels(rsm_ca,fluocats);

Tim = length(rsm_ca);
Ttm = length(accs);
nsynchbins = 10;

where = zeros(length(xcacc)-1,1);


MImat = zeros(nfluocats,naccat);
for k = 1:length(xcacc)-1
    xminus = xcacc(k);
    xplus = xcacc(k+1);
    acclbl = acclabels(xplus-1);
%     
%     tmon = ttime(xminus);
%     tmoff = ttime(xplus);
%    closeon = abs(im4.time-tmon);
%    closeoff = abs(im4.time-tmoff);
     
%     whereon = find(closeon == min(closeon));
%     whereoff = find(closeoff == min(closeoff));
    whereon = xminus;
    whereoff = xplus-1;
%     whereon = whereon(1);
%     whereoff = whereoff(1);
    where(k) = whereon;
    dt = whereoff - whereon + 1;
    
    for i = 1:dt
        imlbl = fluolabels(whereon + i -1);
        MImat(imlbl,acclbl) = MImat(imlbl,acclbl) + 1;
    end
end

MImat = MImat / Tim;

for imlbl = 1:nfluocats
    Pim = length(find(fluolabels == imlbl)) / length(fluolabels);
    for acclbl = 1:naccat
        Ptm = length(find(acclabels == acclbl)) / length(acclabels);
        MImat(imlbl,acclbl) = MImat(imlbl,acclbl)/(Pim*Ptm);
    end
end
clear colormap
figure, colormap('autumn'),imagesc(MImat)
colorbar
ylabel('Fluorescence bins')
xlabel('Acceleration regime')
xticks([1,2,3])
xticklabels({'Deceleration','Ground state','Acceleration'})

%%
rsm_ca = interp1(imtime,cnf.intensity(:,5),ttime);
mm_ca = rsm_ca - min(rsm_ca);
figure, hold on
plot(linspace(1,length(rsm_ca),length(rsm_ca)),mm_ca,'Color',[0 0.1 0.4])
for k = 1:length(where)
    scatter(where(k),mm_ca(where(k)),11,[0 0 0],'filled','MarkerFaceAlpha',0.6), hold on
end
plot(linspace(1,length(rsm_ca),length(accs)),accs*max(mm_ca)+5000,'Color','k','LineWidth',1), hold on

for k = 1:length(xcacc)-1
    xplus = xcacc(k+1);
    xminus = xcacc(k);
    color = acccolors(xplus-1,:);
    
    xl = [xminus, xminus, xplus, xplus];
    yl = [-2000 10000 10000 -2000];
    a = fill(xl,yl,color); hold on,
    a.FaceAlpha = 0.4;
    a.EdgeColor = 'none';
end
axis tight
box off 


%%
clear colormap
nfluocats = 20;
Nimbins = nfluocats;
clear PMI
ttime = tm4.time(1:length(ss_speedM)-1);
imtime = linspace(0,ttime(end),length(m_ca));



for roi = 1:cnf.n_cells
    
    MImat = zeros(Nimbins,Ntmbins);
    rsm_ca = interp1(imtime,cnf.intensity(:,roi),ttime);
    fluocats = define_fluocats(rsm_ca,nfluocats);
    fluolabels = assign_fluolabels(rsm_ca,fluocats);
    
    for k = 1:length(xcacc)-1
        xminus = xcacc(k);
        xplus = xcacc(k+1);
        acclbl = acclabels(xplus-1);

        
        whereon = xminus;
        whereoff = xplus-1;
        dt = whereoff - whereon + 1;

        for i = 1:dt
            imlbl = fluolabels(whereon + i -1);
            MImat(imlbl,acclbl) = MImat(imlbl,acclbl) + 1;
        end
    end

    MImat = MImat / Tim; 
    
    for imlbl = 1:Nimbins
        Pim = length(find(fluolabels == imlbl)) / length(fluolabels);
        for acclbl = 1:Ntmbins
            Ptm = length(find(acclabels == acclbl)) / length(acclabels);
            if Pim>0 && Ptm>0
                MImat(imlbl,acclbl) = MImat(imlbl,acclbl)/(Pim*Ptm);
            else
                MImat(imlbl,acclbl) = 0;
            end
        end
    end
    
        mMI = max(MImat(:));
        MImat = MImat / mMI;
        colormap('gray')
        subplot(10,19,roi), imagesc(MImat), axis off, title(num2str(roi))
        PMI.(['roi_',num2str(roi)]) = MImat;
   
end
suptitle('Cell-resoluted PMI between acceleration and fluorescence')

%% 

PMIcoord = zeros(cnf.n_cells,3);
start = 17 ;
stop = nfluocats;
for roi = 1:cnf.n_cells
    MImat = PMI.(['roi_',num2str(roi)]);
    for j = 1:3
       PMIcoord(roi,j) = mean(MImat(start:stop,j))/mean(MImat(:));
    end
end
figure,
subplot(1,2,1)
scatter3(PMIcoord(:,1),PMIcoord(:,2),PMIcoord(:,3),[],[0 0 0],'filled','MarkerFaceAlpha',0.5), box off
xlabel('Deceleration')
ylabel('Ground state')
zlabel('Acceleration')

K = 3;
indx = kmeans(PMIcoord,K);
cmap = hsv(K);
subplot(1,2,2), hold on
for roi = 1:cnf.n_cells
    ii = indx(roi);
    scatter3(PMIcoord(roi,1),PMIcoord(roi,2),PMIcoord(roi,3),[],cmap(ii,:),'filled','MarkerFaceAlpha',0.5), grid on
    text(PMIcoord(roi,1),PMIcoord(roi,2),PMIcoord(roi,3),num2str(roi),'color',cmap(ii,:),'fontsize',7)
end

xlabel('Deceleration')
ylabel('Ground state')
zlabel('Acceleration')
hold off 

%%
clear colormap
figure, hold on
for k = 1:K
    cl = find(indx == k);
    L = length(cl);
    mintensity = zeros(length(cnf.intensity(:,1)),L);
    mPMI = zeros(nfluocats,3);
    for i = 1:L
        mPMI = mPMI + PMI.(['roi_',num2str(cl(i))]);
        mintensity(:,i) = zero_and_max(cnf.intensity(:,cl(i)).').';  
    end
    mPMI = mPMI/L;
    c = gray(1000);
    ccc = cmap(k,:);
    c = c.*ccc;
    subplot(1,K,k), colormap(gca,c), imagesc(mPMI), box off 
    colorbar 
    ax=gca; 
    ax.YAxis.Color = ccc;
    ax.XAxis.Color = ccc;
    if k == 1
        ylabel('Fluorescence bins')
        xlabel('Acceleration regime')
        xticks([1,2,3])
        xticklabels({'Deceleration','Ground state','Acceleration'})
    else
        xticks([])
        yticks([])
    end
    title(['Cluster ',num2str(k)],'color',ccc)
    clintensity.(['cluster_',num2str(k)]) = mintensity; 
end
hold off 

%%
figure, hold on
S = size(cnf.mask{1,1});
for k = 1:K
    cl = find(indx == k);
    L = length(cl);
    c = cmap(k,:);
    full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
    for i = 1:L
        roi = cl(i);
        posi = cnf.centroid{1,roi};
        x = posi(1);
        y = posi(2);
        %text(x,y,num2str(roi),'Color',c)
        I = cnf.mask{1,roi};
        h = imshow(full); hold on 
        set(h, 'AlphaData', I*0.35) , hold on
    end
end

%%
accounts = [length(find(acclabels==1))/180000,length(find(acclabels==2))/180000,length(find(acclabels==3))/180000]; 

figure, hold on
b = bar([1,2,3],accounts);
b.FaceColor = 'flat';
b.CData(1,:) = [0 0.3 0.9];
b.CData(2,:) = [0.8 0.8 0.8];
b.CData(3,:) = [1 0.4 0.3];
b.FaceAlpha = 0.7;
xticks([])
ylabel('Total time in each state (min)')
xlabel('Acceleration state')


%%
ncat = 1;
speedcats = define_speedcats(ss_speedMs,ncat);
speedlabels = assign_speedlabels(ss_speedMs,speedcats);
xcspd = get_speedregime_boundaries(speedlabels);

plot_all_traces_w_speed(cnf.intensity(:,1:10),xcspd,speedlabels,ss_speedMs,cnf.spikes)

%%
sspikes = zeros(length(cnf.spikes(:,1)),1); 
for k = 1:length(cnf.spikes(:,1))
    sspikes(k) = sum(cnf.spikes(k,:));
end

numcells = [10 20 30 40 50 60 70 80 90];

S = size(cnf.mask{1,1});

ffcells = find(sspikes == 30);  
ffcells = randsample(ffcells,49);
colors = parula(length(numcells));
figure, hold on
oopsc = 1;
for i = 1:length(numcells)
    ffcells = find(sspikes == numcells(i));  
    ffcells = randsample(ffcells,6);
    for k = 1:6
        where = ffcells(k);
        subplot(length(numcells),6,oopsc), hold on
        bkg = cat(3,zeros(S),zeros(S),zeros(S)); %black
        %bkg = cat(3,ones(S),ones(S),ones(S)); %white
        h = imshow(bkg); hold on 
        corrcells = find(cnf.spikes(where,:) == 1); 
        L = length(corrcells);
        c = colors(i,:);
        full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
        for roi = 1:L
          actual = corrcells(roi);
          I = cnf.mask{1,actual};
          h = imshow(full); hold on 
          set(h, 'AlphaData', I*0.9) , hold on
        end
        oopsc = oopsc + 1;
    end
end


%% Subsampling moving bits to check if time spent in each state influences the final result

pairedvals = zeros(2,1); 
c = 1;
for k = 1:length(xcspd)-1
    onset = xcspd(k);
    onlbl = speedlabels(onset); 
    if onlbl == 2
       offset = xcspd(k+1)-1;
       thrshld = (offset - onset + 1)/3000; 
       if thrshld >= 2
           pairedvals(1,c) = onset;
           pairedvals(2,c) = offset;
           c = c + 1;
       end
    end
end

nreps = 36; 
mprobs = zeros(2,nreps);
figure, hold on

for reps = 1:nreps
 
    cumtimemoving = 0;
    possibleidx = linspace(1,length(pairedvals),length(pairedvals)); 
    sampledbits = zeros(2,1);
    c = 1;
    while cumtimemoving < 72 && ~isempty(possibleidx)
        whbit = randsample(possibleidx,1); 
        possibleidx(find(possibleidx == whbit)) = []; 
        onset = pairedvals(1,whbit);
        offset = pairedvals(2,whbit); 
        fullrange = linspace(onset,offset,offset-onset+1); 
        startpt = randsample(fullrange,1); 
        cumtimemoving = cumtimemoving + (offset - startpt + 1)/3000; 
        sampledbits(1,c) = round(startpt/100);
        sampledbits(2,c) = round(offset/100); 
        c = c + 1;
    end

    Pspkloco = zeros(Ncells,1); 
    Pspkrest = zeros(Ncells,1); 
    Nlocobins = 0;
    Nrestbins = 0;

    for k = 1:length(xc)-1
        xminus = xc(k);
        xplus = xc(k+1);
        label = speedlabels(xplus-1);

        tmon = ttime(xminus);
        tmoff = ttime(xplus);
        closeon = abs(imtime-tmon);
        closeoff = abs(imtime-tmoff);

        whereon = round(xminus/100)+1;
        whereoff = round(xplus/100);
        Nb = whereoff-whereon+1;
        if label == 1
            Nrestbins = Nrestbins + Nb;
            for roi = 1:Ncells
                spikes = cnf.spikes(whereon:whereoff,roi);
                whspk = find(spikes==1) + whereon - 1; 
                Nspk = length(whspk);
                Pspkrest(roi) = Pspkrest(roi) + Nspk;
                for i = 1:Nspk
                    coactivated = sum(cnf.spikes(whspk(i),:)) - 1;
                    Nspkrest{roi,:} = cat(2,Nspkrest{roi,:},coactivated); 
                end
            end
        end
    end

    Pspkrest = Pspkrest / Nrestbins;

    Nlocobins = 0;
    for k = 1:length(sampledbits)
        whereon = sampledbits(1,k);
        whereoff = sampledbits(2,k);
        Nb = whereoff - whereon + 1;
        Nlocobins = Nlocobins + Nb;
        for roi = 1:Ncells
            spikes = cnf.spikes(whereon:whereoff,roi);
            whspk = find(spikes == 1) + whereon - 1;
            Nspk = length(whspk);
            Pspkloco(roi) = Pspkloco(roi) + Nspk;
        end
    end

    Pspkloco = Pspkloco / Nlocobins; 
    
    mprobs(1,reps) = mean(Pspkrest);
    mprobs(2,reps) = mean(Pspkloco); 
    
    subplot(6,6,reps)
    meanprobs = [mean(Pspkrest),mean(Pspkloco)];
    incertitude = 2*[std(Pspkrest)/sqrt(Ncells),std(Pspkloco)/sqrt(Ncells)];
    b = bar(meanprobs);hold on
    for roi = 1:Ncells
        plot([1,2],[Pspkrest(roi),Pspkloco(roi)],'Color',[0.8 0.8 0.8 0.5]), hold on
        scatter([1,2],[Pspkrest(roi),Pspkloco(roi)],5,[0.5 0.5 0.5],'filled')
    end
    errorbar(meanprobs,incertitude,'LineStyle','None','Color','k'), 
    b.FaceColor = 'flat';
    b.FaceAlpha = 0.6;
    b.CData(1,:) = colors(1,:);
    b.CData(2,:) = colors(2,:);
    b.EdgeColor = 'none';
    b.BarWidth = 0.9;
    yokpos = [1 7 13 19 25 31];
    xokpos = [31 32 33 34 35 36];
    ylim([0 0.2])
    box off 
    if ~isempty(find(yokpos == reps))
        ylabel('P(spike)')
    else
        yticks([])
    end
    
    if ~isempty(find(xokpos == reps))
        xticks([1 2])
        xticklabels({'R','M'})
    else
        xticks([])
    end
    set(gca,'TickLength',[0,0.01])   
    if reps == nreps
     hfig = gcf;
    end
end
%     figure, hold on
%     meanprobs = [mean(mprobs(1,:)),mean(mprobs(2,:))];
%     incertitude = 2*[std(mprobs(1,:))/sqrt(Ncells),std(mprobs(2,:))/sqrt(Ncells)];
%     LineProps.width = 1;
%     LineProps.edgestyle = '-';
%     LineProps.col = {[0.5 0.5 0.5]};
%     mseb([1,2],meanprobs,incertitude,LineProps,1) 
%     scatter(1,meanprobs(1,1),50,colors(1,:),'filled')
%     scatter(2,meanprobs(1,2),50,colors(2,:),'filled')
%     ylabel('P(spike)')
%     box off 
%     xticklabels({'Resting','Moving'})
%     set(gca,'TickLength',[0,0.01])    

%%

figure, hold on
plot(ss_speedMs*2,'LineWidth',1,'Color','k')
plot(accs-0.2,'LineWidth',1,'Color',[0.2 0.2 0.2])
axis tight
xticks([])
yticks([])
box off




%%
Ncells = cnf.n_cells;
corrMAT = zeros(Ncells,Ncells); 

for k = 1:Ncells
    a = cnf.intensity(:,k);
    for j = 1:Ncells
        b = cnf.intensity(:,j);
        Mc = corrcoef(a,b);
        corrMAT(k,j) = Mc(1,2); 
    end
end

corrMAT = cov(cnf.intensity);
%%
cmap = hot(1000);
%cmap = flipud(cmap);

figure 
reorganised = [];
for j = 1:4
    clnm = activity_clustersK4.clusterregions{:,j}; 
    reorganised = cat(1,reorganised,clnm);
end

orgcorrMAT = zeros(length(reorganised),length(reorganised));
for k = 1:length(reorganised)
    orgcorrMAT(k,:) = corrMAT(reorganised(k),reorganised);
end
colormap(cmap), imagesc(orgcorrMAT), box off 
title('Corrcoefs 4 clusters')



figure
reorganised = [];
for j = 1:3
    clnm = activity_clustersK3.clusterregions{:,j}; 
    reorganised = cat(1,reorganised,clnm);
end

orgcorrMAT = zeros(length(reorganised),length(reorganised));
for k = 1:length(reorganised)
    orgcorrMAT(k,:) = corrMAT(reorganised(k),reorganised);
end
colormap(cmap), imagesc(orgcorrMAT), box off 
title('Corrcoefs 3 clusters')

%% Synchrony in different acceleration states

Nspkcells = zeros(3,cnf.n_cells);
Nbigcoact = zeros(3,cnf.n_cells);
synchregimes = zeros(3,cnf.n_cells); 

figure, hold on
xposi = [2 6 10 14 18 22];
xbars = [1 2 3 ; 5 6 7 ; 9 10 11 ; 13 14 15 ; 17 18 19 ; 21 22 23];
xnames = {'1-10','10-20','20-30','30-40','50-60','>60'};
probchecks = [1 10 ; 10 20 ; 20 30 ;  30 40 ; 50 60 ; 60 170];
[l,~] = size(probchecks); 
for u = 1:l
    synchregimes = zeros(3,cnf.n_cells); 
    for k = 1:length(xcacc)-1
        
        xminus = xcacc(k);
        xplus = xcacc(k+1);
        acclbl = acclabels(xplus-1);

        whereon = round(xminus/100)+1;
        whereoff = round((xplus)/100);

        for j = 1:cnf.n_cells
                cellspk = find(cnf.spikes(whereon:whereoff,j)==1); 
                Nspkcells(acclbl,j) = Nspkcells(acclbl,j) + length(cellspk);
                for cc = 1:length(cellspk)
                    coact = length(find(cnf.spikes(whereon + cellspk(cc) - 1,:) == 1));
                    if coact > probchecks(u,1) && coact < probchecks(u,2)
                        synchregimes(acclbl,j) = synchregimes(acclbl,j) + 1;
                    end
                end
        end

    end

    for k = 1:cnf.n_cells
        for i = 1:3
        synchregimes(i,k) = synchregimes(i,k)/Nspkcells(i,k);
        end
    end
    
    Ncoacc = mean(synchregimes(3,:));
    Ncogr = mean(synchregimes(2,:));
    Ncodec = mean(synchregimes(1,:));
    Nco = [Ncodec, Ncogr, Ncoacc]; 
    Ncoincert = 2*[std(synchregimes(1,:))/sqrt(cnf.n_cells),std(synchregimes(2,:))/sqrt(cnf.n_cells),std(synchregimes(3,:))/sqrt(cnf.n_cells)]; 

    b = bar(xbars(u,:),Nco); hold on
    b.FaceColor = 'flat';
    b.FaceAlpha = 0.7;
    b.CData(3,:) = [0.9216    0.2000    0.1373];
    b.CData(1,:) = [0.4588    1.0157    0.2980];
    b.CData(2,:) = [0    0.1098    0.9608]; 
    errorbar(xbars(u,:),Nco,Ncoincert,'LineStyle','none','Color','k')
end


xticks(xposi)
xticklabels(xnames)
set(gca,'TickLength',[0, 0.05])
xlabel('Number of synchronous dendrites')
ylabel('Probability')
box off

Nspkcells(3,:) = Nspkcells(3,:)*100/length(find(acclabels == 3)); 
Nspkcells(2,:) = Nspkcells(2,:)*100/length(find(acclabels == 2)); 
Nspkcells(1,:) = Nspkcells(1,:)*100/length(find(acclabels == 1)); 
% number of spikes in different conditions
Nspkacc = mean(Nspkcells(3,:));
Nspkgr = mean(Nspkcells(2,:));
Nspkdec = mean(Nspkcells(1,:));
Nspkincert = 2*[std(Nspkcells(1,:))/sqrt(cnf.n_cells), std(Nspkcells(2,:))/sqrt(cnf.n_cells),std(Nspkcells(3,:))/sqrt(cnf.n_cells)];

figure, hold on
b = bar([Nspkdec,Nspkgr,Nspkdec]);
b.FaceColor = 'flat';
b.FaceAlpha = 0.7;
b.CData(3,:) = [0.9216    0.2000    0.1373];
b.CData(1,:) = [0.4588    1.0157    0.2980];
b.CData(2,:) = [0    0.1098    0.9608]; 
errorbar([Nspkdec,Nspkgr,Nspkdec],Nspkincert,'LineStyle','none','Color','k')
xticks([1 2 3])
xticklabels({'Deceleration','Constant','Acceleration'})
xlabel('Regime')
ylabel('Number of spikes')
box off 

% mean number of coactivated cells 

Ncoacc = mean(Nbigcoact(3,:));
Ncogr = mean(Nbigcoact(2,:));
Ncodec = mean(Nbigcoact(1,:));
Nco = [Ncodec, Ncogr, Ncoacc]; 
Ncoincert = 2*[std(Nbigcoact(1,:))/sqrt(cnf.n_cells),std(Nbigcoact(2,:))/sqrt(cnf.n_cells),std(Nbigcoact(3,:))/sqrt(cnf.n_cells)]; 

figure, hold on
b = bar(Nco);
b.FaceColor = 'flat';
b.FaceAlpha = 0.7;
b.CData(3,:) = [0.9216    0.2000    0.1373];
b.CData(1,:) = [0.4588    1.0157    0.2980];
b.CData(2,:) = [0    0.1098    0.9608]; 
errorbar(Nco,Ncoincert,'LineStyle','none','Color','k')
xticks([1 2 3])
xticklabels({'Deceleration','Constant','Acceleration'})
xlabel('Regime')
ylabel('Mean number of coactivated neighbours')
box off 


%%




