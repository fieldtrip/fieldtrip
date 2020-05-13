function [ x ] = ctransform(x)
% CTRANSFORM Copula transformation (empirical CDF)
%   cx = ctransform(x) returns the empirical CDF value along the first
%   axis of x. Data is ranked and scaled within [0 1] (open interval).

[~,x] = sort(x, 1);
[~,x] = sort(x, 1);
x = x / (size(x, 1) + 1);
