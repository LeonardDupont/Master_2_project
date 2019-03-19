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
         for pixel = 1:length(alpha)
            t = alpha(pixel);
            S(1,pixel) = xc + a*cos(t)*cos(tilt) - b*sin(t)*sin(tilt);
            S(2,pixel) = yc + a*cos(t)*sin(tilt) + b*sin(t)*cos(tilt);
         end
    end