function [perimeter] = get_ellipse_perim(geometry)
        
        a = geometry.a;
        b = geometry.b;
        xc = geometry.xc;
        yc = geometry.yc;
        tilt = geometry.tilt;
        
        dangle = 0.001;
        angle = 0; 
        x0 = xc + a*cos(angle)*cos(tilt) - b*sin(angle)*sin(tilt);
        y0 = yc + a*cos(angle)*sin(tilt) + b*sin(angle)*cos(tilt);
        perimeter = 0;
        
        while angle < 2*pi
            angle = angle + dangle;
            x = xc + a*cos(angle)*cos(tilt) - b*sin(angle)*sin(tilt);
            y = yc + a*cos(angle)*sin(tilt) + b*sin(angle)*cos(tilt);
            perimeter = perimeter + sqrt((x0-x)^2 + (y0-y)^2);
            x0 = x;
            y0 = y;
        end
       clear angle 
        
end