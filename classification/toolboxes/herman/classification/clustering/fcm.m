function [C,mu,hist] = fcm(data,numclust,eps,maxiter,exponent,metricfun)

%  FCM implements fuzzy c-means (FCM) clustering algorithm
%  (unsupervised pattern recognition method)
%
%  Used by init_codebook.m or within LVQ object (see lvq.m) OR on its own:
%    [C,mu,hist] = fcm(data,numclust,eps,maxiter,exponent,metricfun)
%
%   INPUT
%       data	 - input features
%       numclust - number of clusters
%       eps      - minimal variation of the optimised objective function that 
%                  determines a stopping criterion (default 1e-6)
%       maxiter  - maximum number of iterations (default 100) 
%       exponent - exponent (commonly 'm') in the objective function (default 2)
%       metricfun- function for calculating the distance metric (between 2 vectors)
%                  (@euclidean_metric is by default) 
%
%   OUTPUT
%       C        - cluster means
%       mu       - membership array [num_datapoints x numclust]
%      hist      - log of the convergence of the objective function 
%
%   SEE ALSO:
%      init_codebook.m
%      lvq.m 
%
%   REFERENCES
%           J. C. Bezdek, Pattern Recognition with Fuzzy Objective Function Algoritms
%           Plenum Press, New York 1981.

% Pawel Herman, 2009

if nargin < 4
    eps = 1e-6;
    maxiter = 100;
end
if nargin < 5
    maxiter = 100;
end
if nargin < 6
    exponent = 2;
end
if nargin < 7
    metricfun = @euclidean_metric;
end

%random initialisation of mu 
mu		= rand(size(data,1),numclust);
rowsum = sum(mu,2);
mu = mu ./ repmat(rowsum,1,numclust);  %normalisation
C		= zeros(numclust,size(data,2));

d = dist(metricfun,data,C);
objfun = sum(sum((d.^2) .* (mu.^exponent)));

niter = 0; 
criterion  = [objfun+10*eps objfun];

while abs(criterion(end-1)-criterion(end)) >= eps && niter < maxiter,
   
   niter = niter + 1;
   
   mf = mu.^exponent;
    
   C = (mf' * data) ./ ( sum(mf)' * ones(1,size(data,2))) ;
       
   d = dist(metricfun,data,C);
   d_exp = d.^(-2/(exponent-1));
   mu = d_exp ./ (sum(d_exp,2) * ones(1,numclust));
    
   objfun = sum(sum((d.^2) .* mf));
   criterion = [criterion objfun];

end

hist = criterion(3:end);


function distance = dist(metricfun,x,y)
% metricfun - metric function defined for two row vectors
% x,y - input vectors (row vectors agreeing on dim) or collection of
%                      vectors - arrays with vectors in rows)
% distance - distance value (scalar) for 2 vectors or array with the
%            distance between corresponding vectors - rows correspond to the first
%            input vector and columns to the second one
distance = zeros(size(x,1),size(y,1));
for i=1:size(x,1)
    for j=1:size(y,1)
        distance(i,j) = feval(metricfun,x(i,:),y(j,:)); 
    end
end

