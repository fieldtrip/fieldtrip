function [y] = randsample(x, k)

% RANDSAMPLE Random sample, with or without replacement. This is a drop-in
% replacement for the matlab funnction with the same name. Not all options
% are supported, an error will be issued if needed.
%
% Y = RANDSAMPLE(N,K) returns Y as a 1-by-K vector of values sampled
% uniformly at random, without replacement, from the integers 1:N.
%
% Y = RANDSAMPLE(POPULATION,K) returns K values sampled uniformly at
% random, without replacement, from the values in the vector POPULATION.
%
% See also RAND, RANDPERM.

% Copyright (C) 2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if nargin>2
  error('only two input variables are supported');
end

if length(x)==1 && isnumeric(x)
  x = 1:x;
end

sel = ceil(rand(1,k)*length(x));
y = x(sel);
