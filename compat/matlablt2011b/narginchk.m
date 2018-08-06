function narginchk(min,max)

% This is a compatibility directory that should only be added to the path on
% MATLAB versions prior to 2011b.
%
% narginchk is not present in older versions.

n = evalin('caller', 'nargin');
assert(n>=min, 'not enough input arguments')
assert(n<=max, 'too many input arguments')

