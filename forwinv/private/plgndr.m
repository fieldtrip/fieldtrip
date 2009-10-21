function [varargout] = funname(varargin)

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
% $Log: plgndr.m,v $
% Revision 1.2  2009/03/12 11:05:03  roboos
% implemented auto-compilation of the mex file in case it is missing
%
% Revision 1.1  2009/01/21 10:32:38  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.3  2008/03/05 16:26:18  roboos
% updated documentation
%
% Revision 1.2  2003/03/11 14:45:37  roberto
% updated help and copyrights
%

% compile the missing mex file on the fly
% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  % try to compile the mex file on the fly
  warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);
  mex(mexsrc);
  cd(pwdir);
  success = true;

catch
  % compilation failed
  disp(lasterr);
  error('could not locate MEX file for %s', mexname);
  cd(pwdir);
  success = false;
end

if success
  % execute the mex file that was juist created
  funname   = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end
