function [varargout] = deepcopy(varargin)

% DEEPCOPY makes a deep copy of an array, and returns a pointer to
% the copy. A deep copy refers to a copy in which all levels of data
% are copied. For example, a deep copy of a cell array copies each
% cell, and the contents of the each cell (if any), and so on.
%
% Example
%   clear a b c
%   a = 1;
%   b = a;            % this is a regular copy
%   c = deepcopy(a);  % this is a deep copy
%   increment(a);     % increment the value of a with one, using pass by reference
%   disp(a);
%   disp(b);
%   disp(c);

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

