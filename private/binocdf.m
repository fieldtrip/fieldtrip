function y=binocdf(x,n,p);

% BINOCDF binomial cumulative distribution function
%
% Y=BINOCDF(X,N,P) returns the binomial cumulative distribution
% function with parameters N and P at the values in X.
%
% See also BINOPDF and STATS (Matlab statistics toolbox)

% compute the cumulative probability for all values up to the maximum
c = cumsum(binopdf(0:max(x(:)),n,p));
y = c(x+1);

% fix rounding errors
y(y<0) = 0;
y(y>1) = 1;

