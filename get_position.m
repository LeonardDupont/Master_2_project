function [x,y] = get_position(angle,geometry)
    
        a = geometry.a;
        b = geometry.b;
        xc = geometry.xc;
        yc = geometry.yc;
        tilt = geometry.tilt;
        
        x = xc + a*cos(angle)*cos(tilt) - b*sin(angle)*sin(tilt);
        y = yc + a*cos(angle)*sin(tilt) + b*sin(angle)*cos(tilt);
        
    
end