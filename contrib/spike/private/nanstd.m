% NANSTD provides a replacement for MATLAB's nanstd that is almost
% compatible.
%
% For usage see STD. Note that the three-argument call with FLAG is not 
% supported.

function Y = nanstd(X, dim)
switch nargin
  case 1
    Y = sqrt(nanvar(X));
  case 2
    Y = sqrt(nanvar(X, dim));
  otherwise
    error('Too many input arguments!');
end
