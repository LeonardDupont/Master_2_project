% Elliptic FISSA - March 2019

clear neuropile
nseg = 5;
neuropile = building_ellipses(cn,'graphics',0,'segment_ellipse',1,'nseg',nseg);

inputfile = '/Users/leonarddupont/Desktop/M2_internship/Code_annex/konkat.tif'; 
neuropile = subR_fluorescence(neuropile,inputfile); 

N = cn.n_cells;

for cell = 42
    
    graphics =1 ;
    if graphics
        cmap = parula(nseg+1);
        figure
        imagesc(neuropile.np_mask_seg{1,cell})
        
        figure, hold on
        
        for seg = 1:nseg+1
            start = 1 + (seg-1)*3 ;
            stop = seg * 3;
            subplot(nseg+1,3,start:stop)
            plot(neuropile.intensity{seg,cell},'color',cmap(seg,:))
            box off
            start = stop + 1;
            stop = stop + 1;
        end    
    end
    
    tic
    disp('1 -- Preparing fluorescence matrix with neuropile subregions --')
    F = initFISSA_mixedF(neuropile,cell);
    toc
    disp('2 -- NNDSVD-based initialisation --')
    [W0,H0] = NNDSVD(F,nseg+1,0); 
    toc

    % . . . . . . . . minimising |F - W*H| according to Frobenius . . . . .
    disp('3 -- Factorising using non-negative constraints --')
    [W,H] = nnmf(F,nseg+1,'w0',W0,'h0',H0); 
    toc
    maxi = max(W(nseg+1,:));
    wc = find(W(nseg+1,:) == maxi);
    dmxd = H(wc,:);
    
    cn.intensity_dm(:,cell) = dmxd; 
    disp(['Cell ',num2str(cell),' done.'])
    
    if graphics 
        figure, hold on
        axh1 = subplot(2,3,1:3);
        mixed = zero_and_max(cn.intensity(:,cell).');
        plot(mixed,'k')
        axis tight, box off
        title('mixed')
        axh2 = subplot(2,3,4:6); 
        demixed = zero_and_max(cn.intensity_dm(:,cell).');
        plot(demixed,'k')
        axis tight, box off 
        title('demixed')
        linkaxes([axh1,axh2],'xy')
        hold off
    end
    
    
    
end





