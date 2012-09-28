function Y = nanstd(X, dim)
switch nargin
  case 1
    Y = sqrt(nanvar(X));
  case 2
    Y = sqrt(nanvar(X, dim));
  otherwise
    error('Too many input arguments!');
end
