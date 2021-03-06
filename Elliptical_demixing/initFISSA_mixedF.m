function [F] = initFISSA_mixedF(neuropile,cell)

nseg = neuropile.nseg;
if neuropile.bkg
    nseg = nseg + 1;
end
t = length(neuropile.intensity{1,cell}); 
F = zeros(nseg+1,t);

for seg = 1:(nseg+1)
    F(seg,:) = neuropile.intensity{seg,cell};
end


end