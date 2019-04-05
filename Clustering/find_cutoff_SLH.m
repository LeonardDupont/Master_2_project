function [Cvopt,Topt,results] = find_cutoff_SLH(Z,distMAT,varargin)
%% April 2019 - CareyLab - leonard.dupont@ens.fr
% .........................................................................
% This function takes the linkage values from the hierarchical clustering
% as well as the distance matrix and finds the optimal cutoff value (here:
% synchrony) according to the Silhouette criterion. This criterion is
% calculated for each point, then averaged over all points of the dataset
% (or of each cluster if this is what one wants). 
%       For each point i of cluster Ci : 
%
%    - a(i) = 1/(Ti -1) * sum( d(i,j) ), j in Ci       
%    will be the overall closeness of point i to other points of Ci. 
%
%    - b(i) = min[i~=j](1/Cj sum( d(i,j) ), j in Cj~=Ci
%    will be the average, maximal closeness of this point to any of the
%    other clusters.
%
%    - s(i) = ( a(i) - b(i) )/ max{a(i),b(i)} is the silhouette of this
%    point. If it indeed belongs to Ci and is 'well clustered', then s(i)
%    will be close to one. Otherwise it will get closer to -1.
%
%  SLH = 1/N * sum(s(i))   where N is the total number of points
%       IS THE FINAL SILHOUETTE INDEX FOR THE CONFIGURATION
%
% .........................................................................
%
% Here, we travel through cutoff values from Cvmax to Cvmin (maximal nb of
% clusters to minimal one) while evaluating the corresponding cluster
% configurations. At the end, the Cv that generated SLH closest to 1 will
% be picked. 
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
parse(ip,varargin{:})
Nmax = ip.Results.Nmax;
Nmin = ip.Results.Nmin;
epsilon = ip.Results.epsilon; 

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
[N,~] = size(distMAT);
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
nth = 3; 
[~,nconfig] = size(clustCONFIG); 
results = zeros(2,nconfig);

for k = 1:nconfig
    
    % first we create a struct that stores the clustering configuration
    T = clustCONFIG(:,k);
    Nclust = clustSIZE(k);
    clear configuration
    skipit = false;
    
    if Nclust == Nmin
        for j=1:Nclust
           cc = length(find(T==j));
           if cc < nth
               results(1,k) = ctf_vals(k);
               results(2,k) = -1e3; 
               skipit = true;
           end
        end
    end
    
    if skipit
        continue
    end
    
    for j = 1:Nclust
        cc = find(T==j);
        configuration.(['c',num2str(j)]) = cc;
    end
    
    SI = zeros(N,1); 
    count = 1;
    for i = 1:Nclust
        
        ci = configuration.(['c',num2str(i)]);
        Ti = length(ci);
        
        
    
        for refpoint = 1:Ti
            
            % - - - - - We start with compactness - - - - - - - - - - - - -
            ai = zeros(Ti-1,1);
            c = 1;
            for ii = 1:Ti
                if ii~=refpoint
                    ai(c) = distMAT(ci(refpoint),ci(ii));
                    c = c+1;
                end
            end
            if Ti > 1
                a = 1/(Ti-1) * sum(ai);
                ai = a ;
                clear a
            else
                ai = 0;
            end
            
            % - - - - - And then closeness to other clusters - - - - - - - 
            howclose = zeros(Nclust-1,1); 
            c = 1;
            for j = 1:Nclust
                if j~=i
                    cj = configuration.(['c',num2str(j)]);
                    Tj = length(cj);
                    close = zeros(Tj,1);
                    for jj = 1:Tj
                        close(jj) = distMAT(cj(jj),ci(refpoint));
                    end
                    
                    howclose(c) = 1/Tj * sum(close); 
                    c = c + 1;
                end
            end
            
            bi = min(howclose); 
            
            if Ti > 1
                si = (bi - ai) / max(ai,bi); 
            else
                si = 0;
            end
            
            SI(count) = si;
            count = count + 1;
        end
    end
    
   SIVal = 1/N * sum(SI); 
   results(1,k) = ctf_vals(k);
   results(2,k) = SIVal;
end


% ----- Now we output the ideal values for the user's comfort! ----------- 
maxS = max(results(2,:)); 
indx = find(results(2,:) == maxS);

Cvopt = results(1,indx);
Topt = clustCONFIG(:,indx); 



end







