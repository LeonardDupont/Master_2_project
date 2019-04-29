function [dFdt] = dF_dt_rechnen(F,time)


dFdt = zeros(length(F)-1,1);
for n = 1:length(F)-1
    dFdt(n) = (F(n+1)-F(n))/(time(n+1)-time(n));
end