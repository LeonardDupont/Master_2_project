function [neuropile] = building_ellipses(cn,varargin)
% March 2019 - Carey lab - leonard.dupont@ens.fr
% .........................................................................
% This function is a first step towards trying to demix signals from 
% adjacent rois created by the Mukamel pipeline. As feared, neighbouring
% units contaminate each other's signal. This programme can be used both to
% implement a suite2P-like approach (unique donut around ellipse) or for a
% FISSA-like approach (donut is divided in nseg of equal area). 
% .........................................................................
% Purkinje cells form elongated ROIs that can be approximated with               
% an ellipse. The programme works as described hereunder:
%
%       1 - Approximates the pixel distribution of the ROI with an
%       ellipse with 95% confidence using covariance (PCA). 
%       2 - Creates a wider ellipse centered on the same origin as the
%       previously-approximated one and finds the pixels that are out 
%       of (1) but within (2). 
%       3 - Generates a mask out of it ---> Suite2P. 
%       4 - Splits the neuropile mask into (nseg) regions ---> FISSA.   
% .........................................................................
%
%  ----- INPUT ----------------------
%
%    cn             carey neuron struct (usual, see documentation)
%                   plz give full struct for the output
%
%    varargin       graphics     @logical    produces plot if 1 (all rois!)
%                   l            @double     additional width of the 2nd
%                                            ellipse
%               segment_ellipse  @logical    1 if you want to segment
%                       nseg     @double     number of segments 
%                       bkg      @logical    if the background should be 
%                                            included
%                       bkgmask  @matrix     cn_bkg.mask basically
%  
%  ----- OUTPUT ---------------------
%
%    neuropile      struct with fields (~similar to cn)   
% .........................................................................


ip = inputParser;
ip.addParameter('graphics',0);
ip.addParameter('l',7);
ip.addParameter('segment_ellipse',0);
ip.addParameter('nseg',5);
ip.addParameter('bkg',0);
ip.addParameter('bkgmask',[]);
parse(ip,varargin{:})
graphics = logical(ip.Results.graphics);
l = ip.Results.l;
segment_ellipse = logical(ip.Results.segment_ellipse);
nseg = ip.Results.nseg;
N = cn.n_cells;
bkg = logical(ip.Results.bkg); 
bkgmask = ip.Results.bkgmask; 

%%  BUILDING THE ELLIPSE AND ITS NEUROPILE DONUT
tic
for roi = 1:N
    
    if rem(roi,25) == 0 || roi == N || roi == 1
        disp(['Building elliptic donut for ROI ',num2str(roi),' out of '...
            ,num2str(N),'.'])
        toc
    end
    
    clear purkinje
    clear geometry
    clear geometry_enlarged
    
    % Calculating eigenvalues, eigen vectors and our new orthogonal base
    themask = cn.mask{1,roi};
    [h,w] = size(themask);
    [xs,ys] = find(themask == 1);  
    xc = mean(xs);
    yc = mean(ys);
    points = [xs ys];
    
    COVM = cov(points);
    [V,D] = eigs(COVM);
    
    d = diag(D);

    geometry.a = sqrt(d(1) * 5.991); % chi-square likelihood at 95% cfdnc
    geometry.b = sqrt(d(2) * 5.991); % p(s < 5.991) = 0.9
    geometry.xc = xc;
    geometry.yc = yc;
    
    %calculating the tilt of the ellipse
    e1 = [1 0];
    if sqrt(d(2)) > sqrt(d(1))
        v1 = V(2,1); %these are x = heights, so actually y in R^2
        v2 = V(2,2); %these are y = width, so actually x in R^2
        geometry.tilt = pi/2 + vector_angles2(e1,[v2,v1]);
        Pb1b2 = [v1 , v2 ; V(1,1) , V(1,2)];
    else
        v1 = V(1,1);
        v2 = V(1,2);
        geometry.tilt = pi/2  + vector_angles2(e1,[v2,v1]);
        Pb1b2 = [v1 , v2 ; V(2,1) , V(2,2)];
    end
 
    
    %this returns the points [x,y] belonging to the perimeter
    Ssmall = border_ellipse(geometry);
 
    
    %geometric properties of the enlarged ellipse are identical, except
    %for semi-axes lengths
    geometry_enlarged.a = geometry.a + l;
    geometry_enlarged.b = geometry.b + l;
    geometry_enlarged.xc = geometry.xc;
    geometry_enlarged.yc = geometry.yc;
    geometry_enlarged.tilt = geometry.tilt;
    
    %set of points for the larger one
    Slarge = border_ellipse(geometry_enlarged);
    
    gradientmask = zeros(h,w);
    np_mask = zeros(h,w);
    
    imx = linspace(1,w,w);
    imy = linspace(1,h,h);
    %creating the mask
    for x=1:length(imx)
        for y=1:length(imy)
            
            insmall = logical(inoutellipse(y,x,geometry) < 1); %within?
            inbig = logical(inoutellipse(y,x,geometry_enlarged) < 1);
            
            if ~insmall && inbig %if outside of small but inside the wider
                np_mask(y,x) = 1;
                gradientmask(y,x) = ...
                                  20*max(0,inoutellipse(y,x,geometry) - 1);
            end
        end
    end
    
 np_mask = logical(np_mask); %to a logical, binary image
 

 maxweight = max(gradientmask(:));
 gradientmask = gradientmask/maxweight;
 neuropile.donutmask{1,roi} = np_mask;
 neuropile.donutgradient{1,roi} = gradientmask;
 
 %% SEGMENTING THE DONUT INTO (nseg) REGIONS 
 
     if segment_ellipse
          Pb2b1 = inv(Pb1b2);
          %this is the transfer matrix from B1 (canonical) to B2 (ellipse)
         % and      [x,y]b2 = Pb2b1 * [x,y]b1' (transfer relation)
         perimeter = get_ellipse_perim(geometry_enlarged);
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
         angle = 0; %rad
         [last_xy(1),last_xy(2)] = get_position(angle,geometry_enlarged);
         angles = zeros(nseg,1);
         angles(end) = 2*pi ;
         
         for s = 2:nseg
             per = 0;
             x0 = last_xy(1);
             y0 = last_xy(2);
             dperim = perimeter/nseg;
             while (per < dperim) && (angle < 2*pi)
                 angle = angle + 0.001;
                 [x,y] = get_position(angle,geometry_enlarged);
                 dis = sqrt((x0-x)^2 + (y0-y)^2);
                 per = per + dis;
                 x0 = x;
                 y0 = y;
             end
             last_xy(1) = x;
             last_xy(2) = y;
             angles(s-1) = angle;
         end
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .    
         np_mask_seg = zeros(h,w);

         for x = 1:h
             for y = 1:w
             % .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
             xy = [x-xc,y-yc];
             xy = Pb2b1 * xy.';
             xy = [xy(2),xy(1)];
             
             theta = vector_angles2([1,0],xy);
             
             
             seg = 1;
             while theta > angles(seg)
                 seg = seg + 1 ; 
             end

            np_mask_seg(x,y) = seg; 
            end
         end
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
     np_mask_seg = np_mask_seg .* np_mask; %convolute to keep the donut only
     themask = (nseg+1)*themask;
     np_mask_seg = np_mask_seg + themask;
     
     if bkg && ~isempty(bkgmask)
         np_mask_seg = np_mask_seg + (nseg+2)*bkgmask;
     end
     
     neuropile.np_mask_seg{1,roi} = np_mask_seg;
     end
     
 %% GRAPHICAL PART
 
    if graphics
        figure;
        ax = gca();
        hold(ax, 'on');
        scatter(ax,ys,xs,'w','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.2),
        hold on
        quiver([yc,yc],[xc,xc],[V(2,1)*sqrt(d(1)),V(2,2)*sqrt(d(2))],...
       [V(1,1)*sqrt(d(1)),V(1,2)*sqrt(d(2))],'color','r','LineWidth',1.5),
       % PUT THIS BACK IF YOU WANT TO PLOT EIGENVECTORS
        hold on
        plot(Ssmall(2,:) , Ssmall(1,:)), hold on
        plot(Slarge(2,:) , Slarge(1,:)), hold on
        if segment_ellipse
            for s = 1:nseg
                angle = angles(s);
                [x,y] = get_position(angle,geometry_enlarged);
                %scatter(y,x,'filled','w'), hold on
            end
            colormap('hot')
            ims = imagesc(ax,np_mask_seg);
            uistack(ims,'bottom')
        end
        imh = image(ax, np_mask);
        %imu = image(ax,gradientmask);
        hold(ax,'off')
        %uistack(imu,'bottom')
        uistack(imh,'bottom')
        xlim([0 w])
        ylim([0,h])
        clear figure 
    end    
end

  neuropile.fov_height = h;
  neuropile.fov_width = w;
  neuropile.n_cells = N;
  neuropile.nseg = nseg;
  if bkg
      neuropile.bkg = true;
  else
      neuropile.bkg = false;
  end

end


    
