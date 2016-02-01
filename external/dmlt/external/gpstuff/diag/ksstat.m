function [snks, snkss] = ksstat(varargin)
%KSSTAT Kolmogorov-Smirnov statistics
%
%   ks = KSSTAT(X) or
%   ks = KSSTAT(X1,X2,...,XJ)
%   returns Kolmogorov-Smirnov statistics in form sqrt(N)*K
%   where M is number of samples.  X is a NxMxJ matrix which
%   contains J MCMC simulations of length N, each with 
%   dimension M. MCMC-simulations can be given as separate
%   arguments X1,X2,... which should have the same length.
%
%   When comparing more than two simulations, maximum of all the
%   comparisons for each dimension is returned (Brooks et al,
%   2003).
%
%   Function returns ks-values in vector KS of length M.
%   An approximation of the 95% quantile for the limiting
%   distribution of sqrt(N)*K with M>=100 is 1.36. ks-values
%   can be compared against this value (Robert & Casella, 2004).
%   In case of comparing several chains, maximum of all the
%   comparisons can be compared to simulated distribution of
%   the maximum of all comparisons obtained using independent 
%   random random numbers (e.g. using randn(size(X))).
%
%   If only one simulation is given, the factor is calculated
%   between first and last third of the chain. 
%
%   Note that for this test samples have to be approximately
%   independent. Use thinning or batching for Markov chains.
%
%   Example: 
%     How to estimate the limiting value when comparing several
%     chains stored in (thinned) variable R. Simulate 100 times
%     independent samples. kss contains then 100 simulations of the
%     maximum of all the comparisons for independent samples.
%     Compare actual value for R to 95%-percentile of kss.
%
%     kss=zeros(1,100);
%     for i1=1:100
%      kss(i1)=ksstat(randn(size(R)));
%     end
%     ks95=prctile(kss,95);
%
%   References:
%    Robert, C. P, and Casella, G. (2004) Monte Carlo Statistical
%      Methods. Springer. p. 468-470.
%    Brooks, S. P., Giudici, P., and Philippe, A. (2003)
%      "Nonparametric Convergence Assessment for MCMC Model
%      Selection". Journal of Computational & Graphical Statistics,
%      12(1):1-22.
%
%   See also
%     PSRF, GEYER_IMSE, THIN

% Copyright (C) 2001-2005 Aki Vehtari
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% In case of one argument split to two halves (first and last thirds)
if nargin==1
  X = varargin{1};
  if size(X,3)==1
    n = floor(size(X,1)/3);
    x = zeros([n size(X,2) 2]);
    x(:,:,1) = X(1:n,:);
    x(:,:,2) = X((end-n+1):end,:);
    X = x;
  end
else
  X = zeros([size(varargin{1}) nargin]);
  for i1=1:nargin
    X(:,:,i1) = varargin{i1};
  end
end

if (size(X,1)<1)
  error('X has zero rows');
end
  
[n1,n2,n3]=size(X);
%if n1<=100
%  warning('Too few samples for reliable analysis');
%end
P = zeros(1,n2);
snkss=zeros(sum(1:(n3-1)),1);
for j1=1:size(X,2)
  ii=0;
  for i1=1:n3-1
    for i2=i1+1:n3
      ii=ii+1;
      snkss(ii,j1)=ksc(X(:,j1,i1),X(:,j1,i2));
    end
  end
  snks(j1)=max(snkss(:,j1));
end

function snks = ksc(x1,x2)
n=numel(x1);
edg=sort([x1; x2]);
c1=histc(x1,edg);
c2=histc(x2,edg);
K=max(abs(cumsum(c1)-cumsum(c2))/n);
snks=sqrt(n)*K;
