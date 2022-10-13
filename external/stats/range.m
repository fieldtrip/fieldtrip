function r = range(x, dim)

% RANGE computes the range (i.e. difference between min and max) for a vector
% or an N-dimensional array. 
%
% Use as
%   r = range(x)
% or you can also specify the dimension along which to look by
%   r = range(x, dim)

% This function was written to be a plugin replacement of the
% similarly-named function in the Matlab stats toolbox.

if nargin < 2
  r = max(x) - min(x);
else
  r = max(x,[],dim) - min(x,[],dim);
end
