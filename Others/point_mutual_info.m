function [PMI] = point_mutual_info(var1,ncats1,var2,ncats2,varargin)
%% April 2019 - CareyLab - leonard.dupont@ens.fr
% .........................................................................
% This function takes in two 1-D vectors specified in var1 - var2 and two
% number of spatial bins ncats1 - ncats2. It then creates ncats categories
% for each variable (respectively) and assigns each datapoint to one of
% them. The programme then outputs a pointwise mutual-information matrix 
% and displays it. If vectors are of different  sizes, the shortest one is
% oversampled to match the length of the other. 
% .........................................................................
%
% - - - - - INPUT - - - - - - - - -
%
%    var1      1-D variable of lenght L1
%    var2      1-D variable of length L2
%    ncats1    number of categories for var1
%    ncats2    number of categories for var2
%  varargin    'withlog'     @logical     PMI = log(PMI) (usual)
%              'plotmat'     @logical     imagesc the matrix
%
% - - - - - OUTPUT - - - - - - - - - 
%
%    PMI       a ncats1 * ncats2 matrix with point-mutual information
%
% .........................................................................
%                               FORMULA 
%
%         PMI(A,B) = p(A,B)/(p(A)p(B)), possibly with a log. 
% 
%   - Here, the probability of any event A, p(A), would be its frequency of
%   occurence in the input data. 
%   - (A,B) is the intersection of events A and B. 
% .........................................................................

ip = inputParser;

ip.addParameter('withlog',1);
ip.addParameter('showmat',1)
parse(ip,varargin{:});

withlog = logical(ip.Results.withlog);
showmat = logical(ip.Results.showmat);

%% If input vectors have different lengths, then we match the size

L1 = length(var1);
L2 = length(var2);

if L1 > L2
    
    q1 = linspace(0,L1,L1);
    q2 = linspace(0,L2,L2);
    var2 = interp1(q2,var2,q1);
    clear q1, clear q2
    
elseif L2 > L1
    
    q1 = linspace(0,L1,L1);
    q2 = linspace(0,L2,L2);
    var1 = interp1(q1,var1,q2);
    clear q1, clear q2
    
end

L = length(var1);
clear L1, clear L2


%% Finding category boundaries for each 1D vector based on (ncats)  
% ------- var1 --------
max1 = max(var1);
min1 = min(var1);
dcat = (max1 - min1)/ncats1; 
cats1 = zeros(ncats1,1);
for k = 1:ncats1
    cats1(k) = k*dcat + min1; 
end

% ------ var2 ---------
max2 = max(var2);
min2 = min(var2);
dcat = (max2 - min2)/ncats2;
cats2 = zeros(ncats2,1);
for j = 1:ncats2
    cats2(j) = j*dcat + min2;
end


%% Assigning labels to each data point of the vectors
% At the end of this part, we have two 1D vectors of size L where label(i)
% is the category ([1:ncat]) to which var(i) belongs. 

labels1 = zeros(L,1);
labels2 = zeros(L,1);

% ------- var1 --------
for k = 1:L
    cat = 1;
    while var1(k) > cats1(cat) && (cat < ncats1)
        cat = cat + 1;
    end
    labels1(k) = cat; 
end
% ------- var2 --------
for j = 1:L
    cat = 1;
    while var2(j) > cats2(cat) && (cat < ncats2)
        cat = cat + 1;
    end
    labels2(j) = cat;
end

%% PMI

PMI = zeros(ncats1,ncats2);

for k = 1:L
    lb1 = labels1(k);
    lb2 = labels2(k);
    PMI(lb1,lb2) = PMI(lb1,lb2) + 1;
end

PMI = PMI / L ; 

for k = 1:ncats1
    P1 = length(find(labels1 == k)) / L;
    for j = 1:ncats2
        P2 = length(find(labels2 == j)) / L;
        if P1>0 && P2>0
            PMI(k,j) = PMI(k,j)/(P1*P2); 
            if withlog
                PMI(k,j) = log(PMI(k,j)); 
            end
        else
            PMI(k,j) = 0;
        end
    end
end

if showmat
    plPMI = PMI / max(PMI(:));
    colormap('parula'), imagesc(plPMI)
    box off 
end


