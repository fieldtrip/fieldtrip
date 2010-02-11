function peerslave(varargin)

% PEERSLAVE starts the low-level peer services and switches to slave mode.
% Subsequently it will wait untill a job comes in and execute it.
%
% Use as
%   peerslave(...)
%
% Optional input arguments should be passed as key-value pairs and can include
%   maxnum     = number (default = inf)
%   maxtime    = number (default = inf)
%   sleep      = number in seconds (default = 0.01)
%   memavail   = number, amount of memory available       (default = inf)
%   cpuavail   = number, speed of the CPU                 (default = inf)
%   timavail   = number, maximum duration of a single job (default = inf)
%   allowhost  = {...}
%   allowuser  = {...}
%   allowgroup = {...}
%   fairshare  = [a, b, c, d]
%   group      = string
%   hostname   = string
%
% See also PEERMASTER, PEERRESET, PEERFEVAL, PEERCELLFUN

% get the optional input arguments
maxnum     = keyval('maxnum',     varargin); if isempty(maxnum),   maxnum=inf; end
maxtime    = keyval('maxtime',    varargin); if isempty(maxtime),  maxtime=inf; end
sleep      = keyval('sleep',      varargin); if isempty(sleep),    sleep=0.01; end
memavail   = keyval('memavail',   varargin);
cpuavail   = keyval('cpuavail',   varargin);
timavail   = keyval('timavail',   varargin);
fairshare  = keyval('fairshare',  varargin);
allowhost  = keyval('allowhost',  varargin); if isempty(allowhost), allowhost = {}; end
allowuser  = keyval('allowuser',  varargin); if isempty(allowuser), allowuser = {}; end
allowgroup = keyval('allowgroup', varargin); if isempty(allowgroup), allowgroup = {}; end
hostname   = keyval('hostname',   varargin);
group      = keyval('group',      varargin);

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
warning off
peer('tcpserver', 'start');
peer('announce',  'start');
peer('discover',  'start');
peer('expire',    'start');
warning on

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

% switch to slave mode
peer('status', 1);

% keep track of the time and number of jobs
stopwatch = tic;
prevtime  = toc(stopwatch);
jobnum    = 0;

while true

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

    try
      % there are many reasons why the execution may fail, hence the elaborate try-catch

      if ~iscell(argin)
        error('input argument should be a cell-array');
      end

      if ~ischar(argin{1})
        error('input argument #1 should be a string');
      end

      fname = argin{1};
      argin = argin(2:end);

      % there are potentially errors to catch from the which() function
      if isempty(which(fname))
        error('Not a valid M-file (%s).', fname);
      end

      % it can be difficult to determine the number of output arguments
      try
        numargout = nargout(fname);
      catch me
        if strcmp(me.identifier, 'MATLAB:narginout:doesNotApply')
          % e.g. in case of nargin('plus')
          numargout = 1;
        else
          rethrow(me);
        end
      end

      if numargout<0
        % the nargout function returns -1 in case of a variable number of output arguments
        numargout = 1;
      end

	  % start measuring the time and memory requirements
	  memprofile on
	  timused = toc(stopwatch);

      % evaluate the function and get the output arguments
	  argout  = cell(1, numargout);
	  [argout{:}] = feval(fname, argin{:});

	  % determine the time and memory requirements
      timused = toc(stopwatch) - timused;
      memstat = memprofile('info');
      memprofile off
      memprofile clear

      % determine the maximum amount of memory that was used during the function evaluation
      memused = max([memstat.mem]) - min([memstat.mem]);

      % Note that the estimated memory is inaccurate, because of
      % the dynamic memory management of Matlab and the garbage
      % collector. Especially on small jobs, the reported memory
      % use does not replect the size of the variables involved in
      % the computation. Matlab is able to squeeze these small jobs
      % in some left-over memory fragment that was not yet deallocated.
      % Larger memory jobs return more reliable measurements.

      fprintf('executing job %d took %f seconds and %d bytes\n', jobnum, timused, memused);

      % collect the output options
      options = {'lastwarn', lastwarn, 'lasterr', lasterr, 'timused', timused, 'memused', memused};

    catch me
      argout  = {};
      % the output options will include the error
      options = {'timused', timused, 'lastwarn', lastwarn, 'lasterr', me};
      % an error was detected while executing the job
      warning('an error was detected during job execution');
    end

    peer('put', joblist.hostid, argout, options, 'jobid', joblist.jobid);
    peer('clear', joblist.jobid);
    clear funname argin argout

  end % isempty(joblist)

end % while true

peer('status', 0);

