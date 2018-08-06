function tf = isequaln(varargin)

% This is a compatibility directory that should only be added to the path on
% MATLAB versions prior to 2012a.
%
% isequalwithequalnans is deprecated in newer MATLAB versions, and
% will be removed, so we strive to only use isequaln in the code.

tf = isequalwithequalnans(varargin{:});

end % function
