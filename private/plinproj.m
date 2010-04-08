function [varargout] = funname(varargin)

% PLINPROJ projects a point onto a line or linepiece
%
% Use as
%   [proj, dist] = plinproj(l1, l2, r, flag)
% where l1 and l2 are the begin and endpoint of the linepiece, and r is 
% the point that is projected onto the line
%
% the optional flag can be:
%   0 (default)  project the point anywhere on the complete line
%   1            project the point within or on the edge of the linepiece

% Copyright (C) 2002-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
    mex -I. -c plinproj.c ; mex plinproj.c plinproj.obj geometry.obj
  else
    mex -I. -c geometry.c
    mex -I. -c plinproj.c ; mex -o plinproj plinproj.o geometry.o
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
