function peerworker(varargin)

% PEERWORKER starts the low-level peer services and switches to worker mode.
% Subsequently it will wait untill a job comes in and execute it.
%
% Use as
%   peerworker(...)
%
% Optional input arguments should be passed as key-value pairs. The
% following options are available to limit the peer network, i.e. to
% form sub-networks.
%   group       = string
%   allowuser   = {...}
%   allowgroup  = {...}
%   allowhost   = {...}
%   refuseuser   = {...}
%   refusegroup  = {...}
%   refusehost   = {...}
% The allow options will prevent peers that do not match the requirements
% to be added to the (dynamic) list of known peers. Consequently, these
% options limit which peers know each other. A controller will not send jobs
% to peers that it does not know. A worker will not accept jobs from a peer
% that it does not know.
%
% The following options are available to limit the number and duration
% of the jobs that the worker will execute.
%   maxnum      = number (default = inf)
%   maxtime     = number (default = inf)
%   maxidle     = number (default = inf)
%
% The following options are available to limit the available resources
% that available for job execution.
%   memavail    = number, amount of memory available       (default = inf)
%   cpuavail    = number, speed of the CPU                 (default = inf)
%   timavail    = number, maximum duration of a single job (default = inf)
%
% See also PEERCONTROLLER, PEERRESET, PEERFEVAL, PEERCELLFUN

% Undocumented options
%   sleep       = number in seconds (default = 0.01)
%   threads     = number, maximum number of threads to use (default = automatic)

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
%
% $Id$
% -----------------------------------------------------------------------

if ft_platform_supports('onCleanup')
  % switch to zombie when finished or when ctrl-c gets pressed
  % the onCleanup function does not exist for older versions
  onCleanup(@peerzombie);
end

% get the optional input arguments
maxnum     = ft_getopt(varargin, 'maxnum', inf);
maxtime    = ft_getopt(varargin, 'maxtime', inf);
maxidle    = ft_getopt(varargin, 'maxidle', inf);
sleep      = ft_getopt(varargin, 'sleep', 0.01);
memavail   = ft_getopt(varargin, 'memavail');
cpuavail   = ft_getopt(varargin, 'cpuavail');
timavail   = ft_getopt(varargin, 'timavail');
threads    = ft_getopt(varargin, 'threads');
group       = ft_getopt(varargin, 'group');
allowuser   = ft_getopt(varargin, 'allowuser', {});
allowgroup  = ft_getopt(varargin, 'allowgroup', {});
allowhost   = ft_getopt(varargin, 'allowhost', {});
refuseuser  = ft_getopt(varargin, 'refuseuser', {});
refusegroup = ft_getopt(varargin, 'refusegroup', {});
refusehost  = ft_getopt(varargin, 'refusehost', {});

if ~isempty(threads) && exist('maxNumCompThreads', 'file')
  % this function is only available from MATLAB version 7.5 (R2007b) upward
  % and has become deprecated in MATLAB version 7.9 (R2009b)
  ws = warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
  maxNumCompThreads(threads);
  warning(ws);
end

% these should be cell arrays
if ~iscell(allowuser) && ischar(allowuser)
  allowuser = {allowuser};
end
if ~iscell(allowgroup) && ischar(allowgroup)
  allowgroup = {allowgroup};
end
if ~iscell(allowhost) && ischar(allowhost)
  allowhost = {allowhost};
end
if ~iscell(refuseuser) && ischar(refuseuser)
  refuseuser = {refuseuser};
end
if ~iscell(refusegroup) && ischar(refusegroup)
  refusegroup = {refusegroup};
end
if ~iscell(refusehost) && ischar(refusehost)
  refusehost = {refusehost};
end

% this should not be user-configurable
% if ~isempty(user)
%   peer('user', user);
% end

% this should not be user-configurable
% if ~isempty(hostname)
%   peer('hostname', hostname);
% end

% the group can be specified by the user
if ~isempty(group)
  peer('group', group);
end

% switch to idle worker mode
peer('status', 2);

% check the current access restrictions
info   = peerinfo;
access = true;
access = access && isequal(info.allowhost,  allowhost);
access = access && isequal(info.allowuser,  allowuser);
access = access && isequal(info.allowgroup, allowgroup);
access = access && isequal(info.refusehost,  refusehost);
access = access && isequal(info.refuseuser,  refuseuser);
access = access && isequal(info.refusegroup, refusegroup);

if ~access
  % impose the updated access restrictions
  peer('allowhost',  allowhost);
  peer('allowuser',  allowuser);
  peer('allowgroup', allowgroup);
  peer('refusehost',  refusehost);
  peer('refuseuser',  refuseuser);
  peer('refusegroup', refusegroup);
end

% check the current status of the maintenance threads
threads = true;
threads = threads && peer('announce', 'status');
threads = threads && peer('discover', 'status');
threads = threads && peer('expire',   'status');
threads = threads && peer('tcpserver', 'status');
% threads = threads && peer('udsserver', 'status');

if ~threads
  % start the maintenance threads
  ws = warning('off');
  peer('announce',  'start');
  peer('discover',  'start');
  peer('expire',    'start');
  peer('tcpserver', 'start');
  % peer('udsserver', 'start');
  warning(ws);
end

% the available resources will be announced and are used to drop requests that are too large
if ~isempty(memavail), peer('memavail', memavail); end
if ~isempty(cpuavail), peer('cpuavail', cpuavail); end
if ~isempty(timavail), peer('timavail', timavail); end

% keep track of the time and number of jobs
stopwatch = tic;
prevtime  = toc(stopwatch);
idlestart = toc(stopwatch);
jobnum    = 0;

while true

  if (toc(stopwatch)-idlestart) >= maxidle
    fprintf('maxidle exceeded, stopping as worker\n');
    break;
  end

  if toc(stopwatch)>=maxtime
    fprintf('maxtime exceeded, stopping as worker\n');
    break;
  end

  if jobnum>=maxnum
    fprintf('maxnum exceeded, stopping as worker\n');
    break;
  end

  joblist = peer('joblist');

  if isempty(joblist)
    % wait a little bit and try again
    pause(sleep);

    % display the time every second
    currtime = toc(stopwatch);
    if (currtime-prevtime>=10)
      prevtime = currtime;
      disp(datestr(now));
    end

  else
    % set the status to "busy worker"
    peer('status', 3);

    % increment the job counter
    jobnum = jobnum + 1;

    % reset the error and warning messages
    lasterr('');
    lastwarn('');

    % get the last job from the list, which will be the oldest
    joblist = joblist(end);

    fprintf('executing job %d from %s@%s (jobid=%d, memreq=%d, timreq=%d)\n', jobnum, joblist.user, joblist.hostname, joblist.jobid, joblist.memreq, joblist.timreq);
    
    % get the input arguments and options
    [argin, options] = peer('get', joblist.jobid);

    % set the options that will be used in the watchdog
    % options = {options{:}, 'controllerid', joblist.hostid}; % add the controllerid as option
    % options = {options{:}, 'timavail', 2*(timavail+1)}; % add the timavail as option, empty is ok

    % evaluate the job
    [argout, options] = peerexec(argin, options);

    % write the results back to the controller
    try
      peer('put', joblist.hostid, argout, options, 'jobid', joblist.jobid);
    catch
      warning('failed to return job results to the controller');
    end

    % remove the job from the tcpserver
    peer('clear', joblist.jobid);

    % remember when the worker becomes idle
    idlestart = toc(stopwatch);

    % set the status to "idle worker"
    peer('status', 2);

  end % isempty(joblist)

end % while true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this masks the regular version, this one only updates the status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function peerzombie
peer('status', 0);

