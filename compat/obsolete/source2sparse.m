function varargout = source2sparse(varargin)

% This function is a backward compatibility wrapper. It allows existing
% MATLAB scripts that do not use the new FieldTrip ft_xxx function naming
% scheme to work with recent versions of the FieldTrip toolbox.
% 
% Please look in ft_xxx for the help of the function that you are looking
% for, where xxx is the name of the function that you were looking for.

% Copyright (C) 2009, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$
prefix    = 'ft_';
funname   = mfilename;
warning([upper(mfilename) ' is only a compatibility wrapper, which will soon be removed. Please instead call ' upper(prefix) upper(funname) '.']);
funhandle = str2func([prefix funname]);
[varargout{1:nargout}] = funhandle(varargin{:});
