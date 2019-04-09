function [] = purkinje_artscape(masks,rois)

if nargin < 2
    N = length(masks);
    rois = linspace(1,N,N);
    cmap = jet(N);
else
    N = length(rois);
    cmap = jet(N);
end

S = size(masks{1,1});
white = cat(3,zeros(S),zeros(S),zeros(S)); 
h = imshow(white); hold on 

for r = 1:N 
   art = rois(r);  
   I = masks{1,art};
   c = cmap(r,:);
   full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
   h = imshow(full); hold on 
   set(h, 'AlphaData', I) , hold on
end

hold off 

end