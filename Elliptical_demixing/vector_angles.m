function [angle] = vector_angles(u,v)
    % given two vectors in the any base, calculates the angle they form.
    % based on u.v = |u| * |v| * cos(u,v). Output in rad. 

        au = abs(u);
        av = abs(v);
        dotuv = dot(au,av);
        normu = norm(u);
        normv = norm(v);

        angle = acos(dotuv/(normu*normv));
        
        v1 = v(1);
        v2 = v(2);
        
        if v2 < 0 && v1 < 0
            angle = pi - angle;
        elseif v2 < 0 && v1 > 0 
            angle = - angle;
        end
end