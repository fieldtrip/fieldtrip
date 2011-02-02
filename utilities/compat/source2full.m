function varargout = funname(varargin)

% This function is a backward compatibility wrapper. It allows existing
% Matlab scripts that do not use the new FieldTrip ft_xxx function naming
% scheme to work with recent versions of the FieldTrip toolbox.
% 
% Please look in ft_xxx for the help of the function that you are looking
% for, where xxx is the name of the function that you were looking for.

% Copyright (C) 2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

prefix    = 'ft_';
funname   = mfilename;
warning([upper(mfilename) ' is only a compatibility wrapper, which will soon be removed. Please instead call ' upper(prefix) upper(funname) '.']);
funhandle = str2func([prefix funname]);
[varargout{1:nargout}] = funhandle(varargin{:});
