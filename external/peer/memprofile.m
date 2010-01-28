function varargout = memprofile(varargin)

%  MEMPROFILE ON starts the profiler and clears previously recorded
%  profile statistics.
% 
%  MEMPROFILE OFF stops the profiler.
% 
%  MEMPROFILE RESUME restarts the profiler without clearing
%  previously recorded memory statistics.
% 
%  MEMPROFILE CLEAR clears all recorded profile statistics.
% 
%  STATS = MEMPROFILE('INFO') suspends the profiler and returns
%  a structure containing the current profiler statistics.

% try to compile the missing mex file on the fly

% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
[p, f, x] = fileparts(funname);
mexsrc = fullfile(p, 'src', [f '.c']);
mexdir = fileparts(funname);

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
  error('could not locate MEX file for %s', mfilename);
  cd(pwdir);
  success = false;
end

if success
  % execute the mex file that was juist created
  funname   = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end

