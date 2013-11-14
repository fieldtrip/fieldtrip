function [s, lam] = covshrinkKPM(x, shrinkvar)
% Shrinkage estimate of a covariance matrix, using optimal shrinkage coefficient.
% INPUT:
% x is n*p data matrix
% shrinkvar : 
% 0: corshrink (default)
% 1: varshrink
%
% OUTPUT:
% s is the posdef p*p cov matrix
% lam is the shrinkage coefficient 
%
% See  J. Schaefer and K. Strimmer.  2005.  A shrinkage approach to 
%   large-scale covariance matrix estimation and implications 
%   for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
% This code is based on their original code http://strimmerlab.org/software.html
% but has been vectorized and simplified by Kevin Murphy.
% Adapted by Marcel van Gerven

if nargin < 2, shrinkvar = 0; end

[n p] = size(x);
if p==1, s=var(x); return; end

switch num2str(shrinkvar)
  
  case '1' % Eqns 10 and 11 of Opgen-Rhein and Strimmer (2007)
    [v, lam] = varshrink(x); 
    dsv = diag(sqrt(v));
    r = corshrink(x);
    s = dsv*r*dsv;
  
  otherwise % case 'D' of Schafer and Strimmer
    v = var(x); 
    dsv = diag(sqrt(v));
    [r, lam] = corshrink(x);
    s = dsv*r*dsv;
    
end
    
%%%%%%%%

function [sv, lambda] = varshrink (x)
% Eqns 10 and 11 of Opgen-Rhein and Strimmer (2007)
[v, vv] = varcov(x);
v = diag(v); vv = diag(vv);
vtarget = median(v);
numerator = sum(vv);
denominator = sum((v-vtarget).^2);
lambda = numerator/denominator;
lambda = min(lambda, 1); lambda = max(lambda, 0);
sv = (1-lambda)*v + lambda*vtarget;
 
function [Rhat, lambda] = corshrink(x)
% Eqn on p4 of Schafer and Strimmer 2005
[n, p] = size(x);
sx = makeMeanZero(x); sx = makeStdOne(sx); % convert S to R
[r, vr] = varcov(sx);
offdiagsumrij2 = sum(sum(tril(r,-1).^2)); 
offdiagsumvrij = sum(sum(tril(vr,-1)));
lambda = offdiagsumvrij/offdiagsumrij2;
lambda = min(lambda, 1); lambda = max(lambda, 0);
Rhat = (1-lambda)*r;
Rhat(logical(eye(p))) = 1;

function [S, VS] = varcov(x)
% s(i,j) = cov X(i,j)
% vs(i,j) = est var s(i,j)
[n,p] = size(x);
xc = makeMeanZero(x); 
S = cov(xc);
XC1 = repmat(reshape(xc', [p 1 n]), [1 p 1]); % size p*p*n !
XC2 = repmat(reshape(xc', [1 p n]),  [p 1 1]); % size p*p*n !
VS = var(XC1 .* XC2, 0,  3) * n/((n-1)^2);

function xc = makeMeanZero(x)
% make column means zero
[n,p] = size(x);
m = mean(x);
xc = x - ones(n, 1)*m; 

function xc = makeStdOne(x)
% make column  variances one
[n,p] = size(x);
sd = ones(n, 1)*std(x);
xc = x ./ sd; 


