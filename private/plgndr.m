function [varargout] = plgndr(varargin)

% PLGNDR associated Legendre function
%
% y = plgndr(n,k,x) computes the values of the associated Legendre functions
% of degree N and order K
%
% implemented as MEX file

% the original implementation was based on "Numerical Recipes in C", version 2.0
% but has been replaced with an equvalent function from GNU Scientific Library

% Copyright (C) 2002, Robert Oostenveld
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

% compile the missing mex file on the fly
% remember the original working directory
pwdir = pwd;
pwdir_ressetter=onCleanup(@()cd(pwdir));

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

mexfullpath=fullfile(mexdir,sprintf('%s.%s',mexname,mexext()));
disp(mexfullpath);
has_mex_func=@()exist(mexfullpath,'file');

if ~has_mex_func()
  % try to compile the mex file on the fly
  ft_warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);
  mex(mexsrc);

  if ~has_mex_func()
    ft_error('could not locate / compile MEX file for %s in %s', ...
                mexname, mexfullpath);
  end
end

% execute the mex file that was just created
funname   = mfilename;
funhandle = str2func(funname);
[varargout{1:nargout}] = funhandle(varargin{:});
