function [Pchance] = est_synchrony_reliability(x,y,taumin)



tx = find(x==1);
ty = find(y==1);
Nx = length(tx);
Ny = length(ty);
N = length(x);

Px = Nx/N;
Py = Ny/N;
Pxy = (2*tau+1)*Px*Py; 

cxy = get_c(tx,ty,taumin);
cyx = get_c(ty,tx,taumin);
 
skss = (cxy + cyx);

Pchance = nchoosek(N,skss) * Pxy^(skss) * (1-Pxy)^(N-skss); 

 function [cxy] = get_c(tx,ty,taumin)
        mx = length(tx);
        my = length(ty);
        cxy = 0;
        for i = 1:mx
            for j = 1:my
                
                dt = tx(i) - ty(j);
                if dt>0 && dt<taumin
                   Jij = 1;
                elseif dt == 0
                   Jij = 1/2;
                else
                   Jij = 0;
                end
               
                cxy = cxy + Jij;
            end
        end
 end


end