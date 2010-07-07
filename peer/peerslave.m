function peerslave(varargin)

% PEERSLAVE starts the low-level peer services and switches to slave mode.
% Subsequently it will wait untill a job comes in and execute it.
%
% Use as
%   peerslave(...)
%
% Optional input arguments should be passed as key-value pairs and can include
%   maxnum      = number (default = inf)
%   maxtime     = number (default = inf)
%   maxidle     = number (default = inf)
%   sleep       = number in seconds (default = 0.01)
%   memavail    = number, amount of memory available       (default = inf)
%   cpuavail    = number, speed of the CPU                 (default = inf)
%   timavail    = number, maximum duration of a single job (default = inf)
%   threads     = number, maximum number of threads to use (default = automatic)
%   allowhost   = {...}
%   allowuser   = {...}
%   allowgroup  = {...}
%   group       = string
%   hostname    = string
%
% See also PEERMASTER, PEERRESET, PEERFEVAL, PEERCELLFUN

% -----------------------------------------------------------------------
% Copyright (C) 2010, Robert Oostenveld
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/
% -----------------------------------------------------------------------

% get the optional input arguments
maxnum     = keyval('maxnum',     varargin); if isempty(maxnum),   maxnum=inf; end
maxtime    = keyval('maxtime',    varargin); if isempty(maxtime),  maxtime=inf; end
maxidle    = keyval('maxidle',    varargin); if isempty(maxidle),  maxidle=inf; end
sleep      = keyval('sleep',      varargin); if isempty(sleep),    sleep=0.01; end
memavail   = keyval('memavail',   varargin);
cpuavail   = keyval('cpuavail',   varargin);
timavail   = keyval('timavail',   varargin);
threads    = keyval('threads',    varargin);
allowhost  = keyval('allowhost',  varargin); if isempty(allowhost), allowhost = {}; end
allowuser  = keyval('allowuser',  varargin); if isempty(allowuser), allowuser = {}; end
allowgroup = keyval('allowgroup', varargin); if isempty(allowgroup), allowgroup = {}; end
group      = keyval('group',      varargin);
hostname   = keyval('hostname',   varargin);

% these should be cell arrays
if ~iscell(allowhost) && ischar(allowhost)
  allowhost = {allowhost};
end
if ~iscell(allowuser) && ischar(allowuser)
  allowuser = {allowuser};
end
if ~iscell(allowgroup) && ischar(allowgroup)
  allowgroup = {allowgroup};
end

% start the maintenance threads
ws = warning('off');
peer('tcpserver', 'start');
peer('announce',  'start');
peer('discover',  'start');
peer('expire',    'start');
warning(ws);

if ~isempty(hostname)
  peer('hostname', hostname);
end

if ~isempty(group)
  peer('group', group);
end

% impose access restrictions
peer('allowhost',  allowhost);
peer('allowuser',  allowuser);
peer('allowgroup', allowgroup);

% the available resources will be announced and are used to drop requests that are too large
if ~isempty(memavail), peer('memavail', memavail); end
if ~isempty(cpuavail), peer('cpuavail', cpuavail); end
if ~isempty(timavail), peer('timavail', timavail); end

if ~isempty(threads) && exist('maxNumCompThreads')
  % this function is only available from Matlab version 7.5 (R2007b) upward
  % and has become deprecated in Matlab version 7.9 (R2009b)
  ws = warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
  maxNumCompThreads(threads);
  warning(ws);
end

% switch to idle slave mode
peer('status', 2);

% keep track of the time and number of jobs
stopwatch = tic;
prevtime  = toc(stopwatch);
idlestart = toc(stopwatch);
jobnum    = 0;

while true

  if (toc(stopwatch)-idlestart) >= maxidle
    fprintf('maxidle exceeded, stopping as slave\n');
    break;
  end

  if toc(stopwatch)>=maxtime
    fprintf('maxtime exceeded, stopping as slave\n');
    break;
  end

  if jobnum>=maxnum
    fprintf('maxnum exceeded, stopping as slave\n');
    break;
  end

  joblist = peer('joblist');

  if isempty(joblist)
    % wait a little bit and try again
    pause(sleep);

    % display the time every second
    currtime = toc(stopwatch);
    if (currtime-prevtime>=1)
      prevtime = currtime;
      disp(datestr(now));
    end

  else
    % set the status to "busy slave"
    peer('status', 3);

    % increment the job counter
    jobnum = jobnum + 1;

    % reset the error and warning messages
    lasterr('');
    lastwarn('');

    % get the last job from the list, which will be the oldest
    joblist = joblist(end);

    fprintf('executing job %d from %s (jobid=%d)\n', jobnum, joblist.hostname, joblist.jobid);
    
    % get the input arguments and options
    [argin, options] = peer('get', joblist.jobid);

    % evaluate the job
    [argout, options] = peerexec(argin, options);

    % write the results back to the master
    try
      peer('put', joblist.hostid, argout, options, 'jobid', joblist.jobid);
    catch
      warning('failed to return job results to the master');
    end

    % remove the job from the tcpserver
    peer('clear', joblist.jobid);

    % remember when the slave becomes idle
    idlestart = toc(stopwatch);

    % set the status to "idle slave"
    peer('status', 2);

  end % isempty(joblist)

end % while true

peer('status', 0);

