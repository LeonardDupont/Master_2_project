function [Q,Qprime,qprime,taumin] = est_trace_synchrony(x,y,varargin)
%% March 2019 - CareyLab - leonard.dupont@ens.fr
% .........................................................................
% This function takes two lists x and y of size t (datapoints acquired at
% a given sampling freq.) with binary values, where 1s occur at times of
% events. These can be peaks or spikes, or anything else you want. This
% algorithm was written fully based on Quiroga et al. (2008) [1] so
% copyright goes to them. 
% .........................................................................
%
%   --- INPUT -----
%
%       x       vector of size t with mx events (1s)
%       y       vector of size t with my events (1s)
%   varargin    'derivative'      @ logical      if you want to compute
%                                                time-resolved synchrony
%
%   --- OUTPUT -----
% 
%        Q     synchrony measure, Q \in [0,1]
%     Qprime   time-resolved synchrony with window size dn
%     qprime   time-resolved antisynchrony (causality)
%     taumin   minimal tau ever defined
% .........................................................................

ip = inputParser;
ip.addParameter('derivative',0);
parse(ip,varargin{:});

derivative = logical(ip.Results.derivative);
%% Q (symmetry)

tx = find(x==1);
ty = find(y==1);

dn = 5;

[cxy,tauminxy] = get_c(tx,ty);
[cyx,tauminyx] = get_c(ty,tx);

taumin = min(tauminxy,tauminyx);

mx = length(tx); %nb of x events
my = length(ty); %nb of y events

if mx~=0 && my~=0
    Q = (cxy + cyx)/sqrt(mx*my);  %simple and straightforward
else
    Q = 0;
end

%% Q'(n)

if derivative
    
    dp = floor(length(x)/dn); %nb of points per window
    aside = rem(length(x),dn); %but this has to be added to the last one
    windows = zeros(dp+1,1);
    windows(1) = 1;
    for w = 1:dp
        windows(w+1) = w * dn + (w==dp)*aside; %stop and start points
    end

    Qval = zeros(dp+1,1);
    qval = zeros(dp+1,1); 
    for w = 1:dp
        cxy = get_c(tx,ty,windows(w+1)); %now we find synchroneous events up
        %to t = windows(w+1). 
        cyx = get_c(ty,tx,windows(w+1));
        Qval(w) = cxy + cyx;
        qval(w) = cxy - cyx;
    end

    Qprime = zeros(dp-1,1);
    qprime = zeros(dp-1,1);
    for w = 1:dp-1

        dQ = Qval(w+1) - Qval(w); %local derivative
        qprime(w) = qval(w+1) - qval(w); 
        
        dnx = length(find(x(windows(w+1):windows(w+2)) == 1)); 
        dny = length(find(y(windows(w+1):windows(w+2)) == 1));
        %finding number of events between the 2 Qvals that we compare

        if dnx*dny ~=0  
            Qprime(w) = dQ/sqrt(dnx*dny); % if dividing is an option
        else
            Qprime(w) = 0;
        end
    end
else
    Qprime = [];
    qprime = []; %if the derivative was not asked for
end

%% Side functions 
% get_c : will calculate the synchrony of the two event trains
% stpfc : simple (actually could use (n>0) to get a logical) function that
%         outputs 1 if x > 0, 0 otherwise. 

    function [cxy,taumin] = get_c(tx,ty,n)
        
        if nargin < 3
            usestep = false;
        else
            usestep = true;
        end
        taumin = 1000;
        mx = length(tx);
        my = length(ty);
        
        cxy = 0;
        
        if mx == 1 && my > 1
            i = 1;
            alli = [];
                for j = 1:my
                    
                    if j == 1
                        allj = [ty(j+1)-ty(j)];
                    elseif j == my 
                        allj = [ty(j)-ty(j-1)];
                    else
                        allj = [ty(j+1)-ty(j),ty(j)-ty(j-1)];
                    end
                    
                    
                    all = cat(2,alli,allj);
                    tau1 = min(all)/2;
                    tau = min(3,tau1);
                    if tau < taumin
                        taumin = tau;
                    end

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

        elseif my == 1 && mx > 1
            j = 1;
            allj = [];
            
            for i = 1:mx
                
                    if i == 1
                        alli = [tx(i+1)-tx(i)];
                    elseif j == my 
                        alli = [tx(i)-tx(i-1)];
                    else
                        alli = [tx(i+1)-tx(i),tx(i)-tx(i-1)];
                    end
                    
                    
                    all = cat(2,alli,allj);
                    tau1 = min(all)/2;
                    tau = min(3,tau1);
                    if tau < taumin
                        taumin = tau;
                    end

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
            
        elseif mx == 1 && my == 1
            dt = tx(1) - ty(1);
            tau = 3;
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
                
            
        elseif mx>1 && my>1
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
                    tau1 = min(all)/2;
                    tau = min(3,tau1);
                    if tau < taumin
                        taumin = tau;
                    end

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
        else %one of the 2 is null
            cxy = 0;
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



