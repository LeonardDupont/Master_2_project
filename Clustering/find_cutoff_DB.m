function [Cvopt,Topt,results] = find_cutoff_DB(Z,distMAT,varargin)
%% April 2019 - CareyLab - leonard.dupont@ens.fr
% .........................................................................
% This function takes the linkage values from the hierarchical clustering
% as well as the distance matrix and finds the optimal cutoff value 
% according to the Davies-Bouldin criterion. This criterion is a
% cluster-dependent value that can be interpreted as:
%
%        DB(cluster) = intra-compactness / inter-dissimilarity
%
% where     intra-compactness : is how spread (around the centroid) the
%           cluster is. 
%           inter-dissimilarity : how far from other clusters it is.
%   
% Actually, the hereabove criterion is calculated for pairs of clusters,
% and for a given cluster the selected value is the highest (worst) one. At
% the end, the sole output value for the cluster configuration is the mean
% of all DBs. The lower, the better. 
% Here, we travel through cutoff values from Cvmax to Cvmin (maximal nb of
% clusters to minimal one) while evaluating the corresponding cluster
% configurations. At the end, the Cv that generated the most compact and
% dissimilar groups wins and is given back to the user.
% .........................................................................
%
%    INPUT
%
%     Z            linkage data
%     distMAT      matrix of distances, symmetrical, real
%     varargin     Nmax          maximal number of clusters (defines Cvmax)
%                  Nmin          minimal number of clusters (defines Cvmin)
%                  epsilon       step used in Cv(n) = Cvmax + n*epsilon
%
%    OUTPUT
%
%     Cvopt         optimal cutoff value (minimises DB criterion)
%     Topt          corresponding tree (actually vector with cl. indices)
%     results       (2*nconfig) vector with all tested Cvopt and
%                   corresponding DB values
%
% 
%  With advice from Jorge Ramirez, PhD. 
% .........................................................................

%% 0 - Initialisation
% Varargin, variables and stuff

ip = inputParser;
ip.addParameter('Nmax',10);
ip.addParameter('Nmin',2);
ip.addParameter('epsilon',1e-4);
ip.addParameter('min_units',1); %default : no restriction 
parse(ip,varargin{:})
Nmax = ip.Results.Nmax;
Nmin = ip.Results.Nmin;
epsilon = ip.Results.epsilon; 
min_units = ip.Results.min_units;

clear cluster
clear configuration
clear ctf_vals
clear clustCONFIG
clear clustSIZE
clear c

clustCONFIG = [];
ctf_vals = [];
clustSIZE = [];
c = 1;

%% 1 - Creating a spectrum of clustering possibilities
% We will travel from Nclust == Nmax to Nclust == Nmin by shifting cutoff
% values (from Cvmax to Cvmin). This will give us a range of clustering 
% configurations which we will then evaluate using the Davies-Bouldin 
% criterion. 

% -------------- Finding the upper bound of the cutoff --------------------
Cv = 0.05;
T = cluster(Z,'Cutoff',Cv,'criterion','distance');
Nclust = max(T(:)); 
while Nclust > Nmax
    Cv = Cv + epsilon;
    T = cluster(Z,'Cutoff',Cv,'criterion','distance');
    Nclust = max(T(:)); 
end
Cvmax = Cv;
ctf_vals(c) = Cvmax;
clustCONFIG(:,c) = T;
clustSIZE(c) = Nclust;

% --- Now create cluster configurations until Nclust = Nmin ---------------
%                  Cv(n) = Cvmax + n*epsilon
c = 2;
while Nclust > Nmin
    Cv = Cv + epsilon;
    T = cluster(Z,'Cutoff',Cv,'criterion','distance');
    nc = max(T(:));
    if nc ~= Nclust %if we reached a lower value of Nclust
        ctf_vals(c) = Cv; %then we store the corresponding cutoff
        clustCONFIG(:,c) = T; %as well as the cluster configuration
        Nclust = nc;
        clustSIZE(c) = Nclust;
        c = c+1;
    end
end
%Cv = Cvmin; just for understandability

%% 2 - Assessing the quality of each clustering configuration
% Now for each of our cutoff values between Cvmax and Cvmin, we will
% evaluate the statistical relevance of our grouping. For this we use
% Davies Bouldin. 

[~,nconfig] = size(clustCONFIG); 
results = zeros(2,nconfig);

for k = 1:nconfig
    
    % first we create a struct that stores the clustering configuration
    T = clustCONFIG(:,k);
    Nclust = clustSIZE(k);
    
    % this whole bit implements the 'min_units' criterion 
    skipit = false;
    badcl = 0;
    for j=1:Nclust
        cc = length(find(T==j));
        if cc < min_units
            badcl = badcl + 1;
        end
    end
    
    Ncheck = Nclust - badcl;
    if Ncheck < Nmin %if Nclusters with more than nth elements is < Nmin
        skipit = true;
    end
    
    if skipit
        results(1,k) = ctf_vals(k);
        results(2,k) = 1e3; 
        continue
    end
    
    %this is the current cluster configuration
    clear configuration
    for j = 1:Nclust
        cc = find(T==j);
        configuration.(['c',num2str(j)]) = cc;
    end
    
    
    % now we can calculate the DB index for this Cv
    Dvals = zeros(1,Nclust); %this will contain the Di index (max(Rij)) for 
    % all clusters (i)
    for i = 1:Nclust
        Rij = zeros(1,Nclust-1); %the list from which we'll pull the max
        cl = 1;
        
        % . . . . . . . . . Spread from the centroid for ci . . . . . . . . 
         ci = configuration.(['c',num2str(i)]);
         Li = length(ci); %number of cells 
         Xi = zeros(Li*(Li-1),1);
         
         count = 1;       
         for ii = 1:Li
             for jj = 1:Li
                if ii ~= jj
                   Xi(count,1) = distMAT(ci(ii),ci(jj));
                   count = count + 1;
                end
             end
         end
                
         Ti = Li*(Li-1);
         Ai = mean(Xi);
         thesum = Xi - Ai;
         thesum = thesum.^2;
         thesum = sum(thesum);

         % ---- actual compactness formula : average distance to Ai -------      
         Si = sqrt(thesum/Ti);
         if isnan(Si)
             Si = 0;
         end
         % ----------------------------------------------------------------
         
         
        for j = 1:Nclust
            if i~=j
                
                % . . . . . . . And now, spread within cj . . . . . . . . .
                cj = configuration.(['c',num2str(j)]);
                Lj = length(cj);
                Xj = zeros(Lj*(Lj-1),1);
                count = 1;
                
                for ii = 1:Lj
                    for jj = 1:Lj
                        if ii ~= jj
                        Xj(count,1) = distMAT(cj(ii),cj(jj));
                        count = count + 1;
                        end
                    end
                end
                
                Tj = Lj*(Lj-1);
                Aj = mean(Xj);
                thesum = Xj - Aj;
                thesum = thesum.^2;
                thesum = sum(thesum);
                
                % -- actual compactness formula : average distance to Aj --
                Sj = sqrt(thesum/Tj);
                if isnan(Sj)
                    Sj = 0;
                end
                % ---------------------------------------------------------
                
                
                % . . . . . Now we deal with INTERcluster distance . . . . 
                Mij = zeros(Li*Lj,1);
                count = 1;
                for ii = 1:Li
                    for jj = 1:Lj
                    Mij(count) = distMAT(ci(ii),cj(jj));
                    count = count+1;
                    end
                end
                
                M = mean(Mij);
                
                % -- Rij = (intra)compactness / (inter)dissimilarity ------
                Rij(cl) = (Si + Sj)/M;
                % ------ Hence, the smaller the better! (in this case) ----
                cl = cl + 1;
            end
        end
        if ~isempty(Rij)
             Dvals(i) = max(Rij); 
        else
            Dvals(i) = 1e3;
        end
        %it's a worst case scenario, we take the greatest Rij among all
        %intercluster comparisons.
    end
    
    DBsum = sum(Dvals);
    DB = (DBsum/Nclust);
    
    results(1,k) = ctf_vals(k);
    results(2,k) = DB;
end


% ----- Now we output the ideal values for the user's comfort! ----------- 
DBs = results(2,:); 
minDB = min(DBs);
indx = find(DBs == minDB);

Cvopt = results(1,indx);
Topt = clustCONFIG(:,indx);


end






