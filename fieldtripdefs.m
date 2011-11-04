function varargout = funname(varargin)

% This function is a backward compatibility wrapper. It allows existing
% Matlab scripts that do not use the new FieldTrip ft_xxx function naming
% scheme to work with recent versions of the FieldTrip toolbox.
% 
% Please look in ft_xxx for the help of the function that you are looking
% for, where xxx is the name of the function that you were looking for.
%
% For this specific function, please look in ft_defaults, not in
% ft_fieldtripdefs or in fieldtripdefs.

% Copyright (C) 2009, Robert Oostenveld
%
% $Id$

funhandle = str2func('ft_defaults');
[varargout{1:nargout}] = funhandle(varargin{:});
