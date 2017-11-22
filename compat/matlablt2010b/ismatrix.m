function status = ismatrix(x)

% This is a compatibility directory that should only be added to the path on
% MATLAB versions prior to 2010b.
%
% ismatrix is not present in older versions.

siz = size(x);
status = numel(siz)==2 && siz(1)>=0 && siz(2)>=0;

end % function


