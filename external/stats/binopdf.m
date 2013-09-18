function y = binopdf(x,n,p)

% BINOPDF binomial probability density function
%
% Y = BINOPDF(X,N,P) returns the binomial probability density
% function with parameters N and P at the values in X.
%
% See also BINOCDF and STATS (Matlab statistics toolbox)

if prod(size(p))>1
  error('probability should be a scalar');
elseif p<0 || p>1
  error('probability should be between 0 and 1');
end

if min(x(:))<0 | max(x(:))>n
  error('X should be in the range 0:N')
end

nk  = gammaln(n + 1) - gammaln(x + 1) - gammaln(n - x + 1);
lny = nk + x.*log( p) + (n - x).*log(1 - p);
y   = exp(lny);

% fix rounding errors
y(y<0) = 0;
y(y>1) = 1;

