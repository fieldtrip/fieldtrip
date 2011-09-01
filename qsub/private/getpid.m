function pid = getpid()
% GETPID returns the process identifier (PID) of the current Matlab
% process.

% this should be implemented in a MEX-file, so throw a warning here
warning('the MEX-implementation of getpid() should be used; returning a surrogate PID');

% and return a surrogate process ID
pid = round(rand(1)*1e7);

end