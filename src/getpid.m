function pid = getpid()

% GETPID returns the process identifier (PID) of the current Matlab
% process.
%
% Use as
%   num = getpid;

% this is to speed up subsequent calls
persistent previous_argout
if ~isempty(previous_argout)
  pid = previous_argout;
  return
end

% this should be implemented in a MEX-file, so throw a warning here
warning('the MEX-implementation of getpid() should be used; returning a surrogate PID');

% and return a surrogate process ID
pid = round(rand(1)*1e7);

% remember for subsequent calls
previous_argout = pid;

