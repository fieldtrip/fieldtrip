function tf = iscolumn(x)

% This is a compatibility directory that should only be added to the path on
% MATLAB versions prior to 2010b.
%
% iscolumn is not present in older versions.

tf = length(size(x))==2 && size(x,2)==1;

end % function
