function nargoutchk(min,max)

% This is a compatibility directory that should only be added to the path on
% MATLAB versions prior to 2011b.
%
% nargoutchk is not present in older versions.

n = evalin('caller', 'nargout')
assert(n>=min, 'not enough output arguments')
assert(n<=max, 'too many output arguments')

