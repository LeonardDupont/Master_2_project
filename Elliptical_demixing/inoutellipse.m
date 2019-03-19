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