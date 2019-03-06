function [neuropile] = building_ellipses(cn,varargin)
% March 2019 - Carey lab - leonard.dupont@ens.fr
% .........................................................................
% This function is a first step towards trying to demix signals from 
% adjacent rois created by the Mukamel pipeline. As feared, neighbouring
% units contaminate each other's signal. Here, we use a suite2P - like 
% strategy to generate a neuropile mask around each roi which we will 
% further use to calculate a mean, surrounding fluorescence for each cell
% in each frame. This shall hopefully demix the signals.
% .........................................................................
% Purkinje cells form elongated ROIs that can be approximated with               
% an ellipse. The programme works as described hereunder:
%
%       1 - Approximates the pixel distribution of the ROI with an
%       ellipse with 95% confidence using covariance (PCA). 
%       2 - Creates a wider ellipse centered on the same origin as the
%       previously-approximated one and finds the pixels that are out 
%       of (1) but within (2). 
%       3 - Generates a mask out of it.
%       4 - (optional) Producing spatial weighing of the neuropile. 
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
%  
%  ----- OUTPUT ---------------------
%
%    neuropile      struct with fields (similar to cn)   
% .........................................................................


ip = inputParser;
ip.addParameter('graphics',0);
ip.addParameter('l',7)
parse(ip,varargin{:})
graphics = logical(ip.Results.graphics);
l = ip.Results.l;

N = cn.n_cells;

tic
for roi = 1:N
    
    if rem(roi,25) == 0 | roi == N
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
    
    %calculating ellipse rotation angle
    e1 = [1 0]; %angle 0 is horizontal
    if sqrt(d(2)) > sqrt(d(1))
        geometry.tilt = vector_angles(e1,V(2,:));
    else
        geometry.tilt = vector_angles(e1,V(1,:));
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
 neuropile.fov_height = h;
 neuropile.fov_width = w;
 neuropile.n_cells = N;
 neuropile.donutmask{1,roi} = np_mask;
 neuropile.donutgradient{1,roi} = gradientmask;
 
 
 
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
        imh = image(ax, np_mask);
        imu = image(ax,gradientmask);
        hold(ax,'off')
        uistack(imu,'bottom')
        uistack(imh,'bottom')
        xlim([0 w])
        ylim([0,h])
        clear figure 
    end
    

    
end


%% annex functions

    function [S] = border_ellipse(geometry)
    %given geometry properties of the ellipse, 
    %gives back the list of coordinates for the 
    %points forming its perimeter. 
    
    
         a = geometry.a;
         b = geometry.b;
         xc = geometry.xc;
         yc = geometry.yc;
         tilt = geometry.tilt;

         alpha = linspace(0,2*pi,2*pi/0.01);

         S = zeros(2,length(alpha));
         for k = 1:length(alpha)
            t = alpha(k);
            S(1,k) = xc + a*cos(t)*cos(tilt) - b*sin(t)*sin(tilt);
            S(2,k) = yc + a*cos(t)*sin(tilt) + b*sin(t)*cos(tilt);
         end
    end



    function [inout] = inoutellipse(x,y,geometry)
    % using the cartesian equation of the ellipse, returns a value
    % that will be < 1 if the point is inside of it (1 if on the perimeter)
    % and > 1 otherwise.
   
    
        a = geometry.a;
        b = geometry.b;
        xc = geometry.xc;
        yc = geometry.yc;
        tilt = geometry.tilt;


        inout = ( ((x-xc)*cos(tilt) + (y-yc)*sin(tilt))/a)^2 + ...
        ( ((x-xc)*sin(tilt) - (y-yc)*cos(tilt))/b)^2;
    end

    function [angle] = vector_angles(u,v)
    % given two vectors in the any base, calculates the angle they form.
    % based on u?v = |u| * |v| * cos(u,v). Output in rad. 

        dotuv = dot(u,v);
        normu = norm(u);
        normv = norm(v);

        angle = acos(dotuv/(normu*normv));
        
        v1 = v(1);
        v2 = v(2);
        
        if (v1 > 0) && (v2 < 0)
            angle = 2*pi - angle;
        elseif (v1 < 0) && (v2 > 0)
            angle = pi - angle;
        elseif (v1 < 0) && (v2 < 0)
            angle = pi + angle;
        end  
    end


end