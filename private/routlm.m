function [varargout] = routlm(varargin)

% ROUTLM computes the projection of a point from its la/mu parameters
% these equal the "Barycentric" coordinates
%
% Use as
%   [proj] = routlm(v1, v2, v3, la, mu)
% where v1, v2 and v3 are three vertices of the triangle

% Copyright (C) 2002-2009, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the native MATLAB implementation.
% The mex file is many times faster and is therefore preferred.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [proj] = routlm(v1, v2, v3, la, mu);
% proj = (1-la-mu)*v1 + la*v2 + mu*v3;

% compile the missing mex file on the fly
% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  % try to compile the mex file on the fly
  ft_warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);
  
  if ispc
    mex -I. -c geometry.c
    mex -I. -c routlm.c ; mex routlm.c routlm.obj geometry.obj
  else
    mex -I. -c geometry.c
    mex -I. -c routlm.c ; mex -o routlm routlm.o geometry.o
  end
  
  cd(pwdir);
  success = true;

catch
  % compilation failed
  disp(lasterr);
  ft_error('could not locate MEX file for %s', mexname);
  cd(pwdir);
  success = false;
end

if success
  % execute the mex file that was juist created
  funname   = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end
