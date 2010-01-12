function peerslave(varargin)

% PEERSLAVE starts the low-level peer services and switches to slave mode.
% Subsequently it will wait untill a job comes in and execute it.
%
% Use as
%   peerslave(...)
%
% Optional input arguments should be passed as key-value pairs and can include
%   maxnum   = number (default = inf)
%   maxtime  = number (default = inf)
%   sleep    = number in seconds (default = 0.01)
%
% See also PEERMASTER, PEERRESET, PEERFEVAL, PEERCELLFUN

% get the optional input arguments
maxnum  = keyval('maxnum',  varargin); if isempty(maxnum), maxnum=inf; end
maxtime = keyval('maxtime', varargin); if isempty(maxtime), maxtime=inf; end
sleep   = keyval('sleep',   varargin); if isempty(sleep), sleep=0.01; end

warning off
peer('tcpserver', 'start');
peer('announce',  'start');
peer('discover',  'start');
peer('expire',    'start');
warning on

peer('status', 1);

% keep track of the time and number of jobs
stopwatch = tic;
prevtime  = toc(stopwatch);
jobnum = 0;

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

      argout  = cell(1, numargout);
      elapsed = toc(stopwatch);
      [argout{:}] = feval(fname, argin{:});
      elapsed = toc(stopwatch) - elapsed;
      fprintf('executing job %d took %f seconds\n', jobnum, elapsed);

      % collect the output options
      options = {'elapsed', elapsed, 'lastwarn', lastwarn, 'lasterr', lasterr};

    catch
      argout  = {};
      % the output options will include the error
      options = {'elapsed', elapsed, 'lastwarn', lastwarn, 'lasterr', lasterr};
      % an error was detected while executing the job
      warning('an error was detected during job execution');
    end

    peer('put', joblist.hostid, argout, options, joblist.jobid);
    peer('clear', joblist.jobid);
    clear funname argin argout

  end % isempty(joblist)

end % while true

peer('status', 0);

