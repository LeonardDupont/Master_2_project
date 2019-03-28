function [Q,Qprime,qprime] = est_trace_synchrony(x,y,varargin)
%% March 2019 - CareyLab - leonard.dupont@ens.fr
% .........................................................................
% This function takes two fluorescence traces {t1,t2} and evaluates their
% synchrony using a sliding, rectangular window of size (window) seconds.
% The formulas are taken from Quiroga et al. (2008) and the algorithm
% adapted to our purpose.
% .........................................................................
%
%   --- INPUT -----
%
%       x
%       y
%
%
%   --- OUTPUT -----
% 
%              
% .........................................................................

ip = inputParser;
ip.addParameter('derivative',0);
parse(ip,varargin{:});

derivative = logical(ip.Results.derivative);
%% Q (symmetry)

tx = find(x==1);
ty = find(y==1);
mx = length(tx);
my = length(ty);

dn = 5;

cxy = get_c(tx,ty);
cyx = get_c(ty,tx);

Q = (cxy + cyx)/sqrt(mx*my); 

%% Q'(n)

if derivative
    
    dp = floor(length(x)/dn);
    aside = rem(length(x),dn);
    windows = zeros(dp+1,1);
    windows(1) = 1;
    for w = 1:dp
        windows(w+1) = w * dn + (w==dp)*aside;
    end

    Qval = zeros(dp+1,1);
    qval = zeros(dp+1,1); 
    for w = 1:dp
        cxy = get_c(tx,ty,windows(w+1));
        cyx = get_c(ty,tx,windows(w+1));
        Qval(w) = cxy + cyx;
        qval(w) = cxy - cyx;
    end

    Qprime = zeros(dp-1,1);
    qprime = zeros(dp-1,1);
    for w = 1:dp-1

        dQ = Qval(w+1) - Qval(w);
        qprime(w) = qval(w+1) - qval(w); 
        dnx = length(find(x(windows(w+1):windows(w+2)) == 1)); 
        dny = length(find(y(windows(w+1):windows(w+2)) == 1));

        if dnx*dny ~=0  
            Qprime(w) = dQ/sqrt(dnx*dny);  
        else
            Qprime(w) = 0;
        end
    end
else
    Qprime = [];
    qprime = [];
end



    function [cxy] = get_c(tx,ty,n)
        
        if nargin < 3
            usestep = false;
        else
            usestep = true;
        end
        
        mx = length(tx);
        my = length(ty);
        cxy = 0;
        for i = 1:mx
             if i == 1
                    alli = [tx(i+1)-tx(i)];
                elseif i == mx
                    alli = [tx(i)-tx(i-1)];
                else
                    alli = [tx(i+1)-tx(i),tx(i)-tx(i-1)];
             end
            for j = 1:my
                
                if j == 1
                    allj = [ty(j+1)-ty(j)];
                elseif j == my
                    allj = [ty(j)-ty(j-1)];
                else
                    allj = [ty(j+1)-ty(j),ty(j)-ty(j-1)];
                end
                
                all = cat(2,alli,allj);
                tau = min(all)/2;
                
                dt = tx(i) - ty(j);
                if dt>0 && dt<tau
                   Jij = 1;
                elseif dt == 0
                   Jij = 1/2;
                else
                   Jij = 0;
                end
                
                if usestep
                    step = stpfc(n-tx(i));
                else
                    step = 1;
                end
                
                cxy = cxy + Jij*step;
            end
        end
    end
   
    function [step] = stpfc(x)
        if x > 0
            step = 1;
        else
            step = 0;
        end
    end


end



