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

if nargin == 2
   others = linspace(1,length(masks),length(masks));
   others(rois) = [];
   for roth = 1:length(others)
       art = others(roth);
       I = masks{1,art};
       c = [1 1 1];
       full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3));
       h = imshow(full); hold on
       set(h,'AlphaData',I*0.1)
   end
end

hold off 

end