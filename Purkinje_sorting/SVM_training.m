%% March 20th, 2019 - Carey lab - leonard.dupont@ens.fr
%  Supersiving users : Diogo, Hugo, Jorge, Leonard
%  Building an SVM classifier based on mask shape for Purkinje refinement
%  post Mukamel segmentation. 
%% Extracting shape criterion
M = cn2.n_cells;
N = cn.n_cells;

criteria1 = extract_mask_criteria(cn);
criteria2 = extract_mask_criteria(cn2);

deviations = cat(1,criteria1(:,1),criteria2(:,1)); 
spreads = cat(1,criteria1(:,2),criteria2(:,2));
elongations = cat(1,criteria1(:,3),criteria2(:,3));
smoothness = cat(1,criteria1(:,4),criteria2(:,4)); 


bigcn.n_cells = M+N;
for roi = 1:N
    bigcn.mask{1,roi} = cn.mask{1,roi};
end
for roi = 1:M
    bigcn.mask{1,N+roi} = cn2.mask{1,roi};
end

%% Building supervised-learning dataset
N = bigcn.n_cells;
users = {'Leo','Hugo','Jorge','Diogo'};
nuser = length(users); 
good_purkinje_all = zeros(nuser,N);
user = 1; 

clear good_purkinje
% 1 - manual sorting
manual_roi_sorting(bigcn);

% 2 - get the sorting vector
good_purkinje = getglobal_purkinje;
good_purkinje_all(user,:) = good_purkinje;
user = user + 1;


good_pkj_stat = mean(good_purkinje_all,1); 
good_pkj_stat(good_pkj_stat < 0.6) = 0;
good_pkj_stat(good_pkj_stat >= 0.6) = 1;



 
%%

response = good_pkj_stat;
T = ...
  table(deviations,spreads,elongations,smoothness,response.', ...
  'VariableNames',{'deviations','spreads','elongations','smoothness',...
  'response'}); 
clear good_pkj_stat
clear good_purkinje_all
clear good_purkinje

classificationLearner

%%
good_deviations = deviations(good_pkj_stat==1);
bad_deviations = deviations(good_pkj_stat==0);

good_elongations = elongations(good_pkj_stat==1);
bad_elongations = elongations(good_pkj_stat==0);

good_spreads = spreads(good_pkj_stat==1);
bad_spreads = spreads(good_pkj_stat==0);

good_smoothness = smoothness(good_pkj_stat==1);
bad_smoothness = smoothness(good_pkj_stat==0); 



figure, hold on
scatter3(good_deviations,good_elongations,good_spreads,[],[0.6 0.2 0.8],'filled','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',1), hold on
scatter3(bad_deviations,bad_elongations,bad_spreads,[],[1 0.2 0.3],'filled','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',1)
xlabel('deviation')
ylabel('elongation')
zlabel('spread')
legend('Real cells','Fake cells')
grid on, hold off 

%% 
save('SVM_Pkj.mat','Pkj_sorter')

