function [wfreqs] = slidingwindow_freq(binaryspikes)

T = length(binaryspikes);
fs = 30; 
wsize = 0.5; %s
Twsize = round(wsize * fs);

Nit = floor(T/Twsize);
aside = rem(T,Twsize);
wlimits = zeros(Nit,1); 
for k = 1:Nit
    wlimits(k) = (k-1)*Twsize + 1; 
end
wlimits(end) = wlimits(end) + aside;

wfreqs = zeros(Nit-1,1); 

for k = 1:Nit-1
    start = wlimits(k);
    stop = wlimits(k+1);
    Nspk = length(find(binaryspikes(start:stop) == 1));
    wfreqs(k) = Nspk*fs/(stop-start); 
end