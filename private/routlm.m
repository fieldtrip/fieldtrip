function [varargout] = funname(varargin)

% ROUTLM computes the projection of a point from its la/mu parameters
% these equal the "Barycentric" coordinates
%
% Use as
%   [proj] = routlm(v1, v2, v3, la, mu)
% where v1, v2 and v3 are three vertices of the triangle

% Copyright (C) 2002-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the native Matlab implementation.
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
  warning('trying to compile MEX file from %s', mexsrc);
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
