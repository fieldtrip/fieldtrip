function mutexlock(lockfile, timeout)

% MUTEXLOCK creates a lockfile, or if it already exists, waits until
% another process removes the lockfile and then creates it. This function
% can be used for "mutual exclusion", i.e. executing multiple processes in
% parallel where part of the processing is not allowed to run
% simultaneously.
%
% Use as
%   mutexlock(lockfile, timeout)
%
% See also MUTEXUNLOCK and http://en.wikipedia.org/wiki/Mutual_exclusion

if isempty(lockfile)
  % no lockfile was specified, return immediately
  return
end

if nargin<2
  timeout = 3600; % seconds
end

stopwatch = tic;
while exist(lockfile, 'file')
  % wait untill another process removes the lockfile
  pause(1);
  if toc(stopwatch)>timeout
    error('timeout exceeded waiting for the lockfile to be removed');
  end
end

% create the lockfile
fid = fopen(lockfile, 'wb');
if fid<0
  error('cannot open lockfile');
end
flose(fid);
