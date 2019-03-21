function [] = purkinje_artscape(masks)

N = length(masks);
S = size(masks{1,1});
cmap = parula(N);


white = cat(3,zeros(S),zeros(S),zeros(S)); 
h = imshow(white); hold on 

for cell = 1:N 
   I = masks{1,cell};
   c = cmap(cell,:);
   full = cat(3,ones(S)*c(1),ones(S)*c(2),ones(S)*c(3)); 
   h = imshow(full); hold on 
   set(h, 'AlphaData', I) , hold on
end

hold off 

end