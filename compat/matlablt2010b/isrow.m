function tf = isrow(x)

% This is a compatibility directory that should only be added to the path on
% MATLAB versions prior to 2010b.
%
% isrow is not present in older version.

tf = length(size(x))==2 && size(x,1)==1;

end % function
