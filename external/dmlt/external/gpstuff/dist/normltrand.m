function result = normltrand(mu,sigma2,left)
%NORMLTRAND  compute random draws from a left-truncated normal
%            distribution, with mean = mu, variance = sigma2
% ------------------------------------------------------
% USAGE: y = normltrand(mu,sigma2,left)
% where:   mu = mean (scalar or vector)
%      sigma2 = variance (scalar or vector)
%        left = left truncation point (scalar or vector)
% ------------------------------------------------------
% RETURNS: y = (scalar or vector) the size of mu, sigma2
% ------------------------------------------------------
% NOTES: This is merely a convenience function that
%        calls normtrand with the appropriate arguments
% ------------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

% Anyone is free to use these routines, no attribution (or blame)
% need be placed on the author/authors.

% 2009-01-08 Aki Vehtari - Fixed Naming

if nargin ~= 3
error('normltrand: Wrong # of input arguments');
end;

right = mu + 5*sqrt(sigma2);

result = normtrand(mu,sigma2,left,right);
