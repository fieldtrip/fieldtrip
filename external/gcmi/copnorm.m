function x = copnorm(x)
% COPNORM Copula normalisation
%   cx = copnorm(x) returns standard normal samples with the same empirical
%   CDF value as the input. Operates along the first axis.
%   Equivalent to cx = norminv(ctransform(x))
%
[~,x] = sort(x, 1);
[~,x] = sort(x, 1);
x = x / (size(x, 1) + 1);
x = -sqrt(2).*erfcinv(2*x);
