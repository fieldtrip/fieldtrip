function varargout = peercellfun(fname, varargin)

% PEERCELLFUN apply a function to each element of a cell-array. The
% function execution is done in parallel on all avaialble peers.
%
% Use as
%   argout = peercellfun(fname, x1, x2, ...)
%
% This function has a number of optional arguments that have to passed
% as key-value pairs at the end of the list of input arguments. All other
% input arguments (including other key-value pairs) will be passed to the
% function to be evaluated.
%   UniformOutput  = boolean (default = false)
%   StopOnError    = boolean (default = true)
%   ResubmitTime   = number, amount of time for a job to be resubmitted (default is automatic)
%   MaxBusy        = number, amount of slaves allowed to be busy (default = inf)
%   diary          = string, can be 'always', 'never', 'warning', 'error' (default = 'error')
%   memreq         = number
%   timreq         = number
%   sleep          = number
%   order          = string, can be 'random' or 'original' (default = 'random')
%
% Example
%   fname = 'power';
%   x1    = {1, 2, 3, 4, 5};
%   x2    = {2, 2, 2, 2, 2};
%   y     = peercellfun(fname, x1, x2);
%
% See also PEERMASTER, PEERSLAVE, PEERLIST, PEERINFO, PEERFEVAL, CELLFUN, DFEVAL

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

if matlabversion>=7.8
  % switch to zombie when finished or when ctrl-c gets pressed
  % the onCleanup function does not exist for older versions
  onCleanup(@peerzombie);
end

% locate the begin of the optional key-value arguments
optbeg = find(cellfun(@ischar, varargin));
optarg = varargin(optbeg:end);

% assert that the user does not specify obsolete options
keyvalcheck(optarg, 'forbidden', {'timcv'});

% get the optional input arguments
UniformOutput = keyval('UniformOutput', optarg); if isempty(UniformOutput), UniformOutput=false; end
StopOnError   = keyval('StopOnError',   optarg); if isempty(StopOnError),   StopOnError=true;    end
ResubmitTime  = keyval('ResubmitTime',  optarg);
MaxBusy       = keyval('MaxBusy',       optarg); if isempty(MaxBusy),       MaxBusy=inf;         end
memreq        = keyval('memreq',        optarg); if isempty(memreq),        memreq=1024^3;       end % assume 1 GB
timreq        = keyval('timreq',        optarg); if isempty(timreq),        timreq=3600;         end % assume 1 hour
sleep         = keyval('sleep',         optarg); if isempty(sleep),         sleep=0.05;          end
diary         = keyval('diary',         optarg); if isempty(diary),         diary='error';       end % 'always', 'never', 'warning', 'error' 
order         = keyval('order',         optarg); if isempty(order),         order='random';      end % 'random', 'original'

% convert from 'yes'/'no' into boolean value
UniformOutput = istrue(UniformOutput);

% skip the optional key-value arguments
if ~isempty(optbeg)
  varargin = varargin(1:(optbeg-1));
end

if isa(fname, 'function_handle')
  % convert the function handle back into a string (e.g. @plus should be 'plus')
  fname = func2str(fname);
end

% there are potentially errors to catch from the which() function
if isempty(which(fname))
  error('Not a valid M-file (%s).', fname);
end

% determine the number of input arguments and the number of jobs
numargin    = numel(varargin);
numjob      = numel(varargin{1});

% it can be difficult to determine the number of output arguments
try
  numargout = nargout(fname);
catch nargout_err
  if strcmp(nargout_err.identifier, 'MATLAB:narginout:doesNotApply')
    % e.g. in case of nargin('plus')
    numargout = 1;
  else
    rethrow(nargout_err);
  end
end

if numargout<0
  % the nargout function returns -1 in case of a variable number of output arguments
  numargout = 1;
end

% check the input arguments
for i=1:numargin
  if ~isa(varargin{i}, 'cell')
    error('input argument #%d shoudl be a cell-array', i+1);
  end
  if numel(varargin{i})~=numjob
    error('inconsistent number of elements in input #%d', i+1);
  end
end

% check the availability of peer slaves
list = peerlist;
list = list([list.status]==2 | [list.status]==3);
if isempty(list)
  warning('there is no peer available as slave, reverting to local cellfun');
  % prepare the output arguments
  varargout = cell(1,numargout);
  % use the standard cellfun
  [varargout{:}] = cellfun(str2func(fname), varargin{:}, 'UniformOutput', UniformOutput);
  return
end

% prepare some arrays that are used for bookkeeping
jobid       = nan(1, numjob);
puttime     = nan(1, numjob);
timused     = nan(1, numjob);
memused     = nan(1, numjob);
submitted   = false(1, numjob);
collected   = false(1, numjob);
busy        = false(1, numjob);
submittime  = inf(1, numjob);
collecttime = inf(1, numjob);
resubmitted = [];    % this will contain a growing list with structures

% remove any remains from an aborted previous call
joblist = peer('joblist');
for i=1:length(joblist)
  peer('clear', joblist(i).jobid);
end

% start the timer
stopwatch = tic;

% these are used for printing feedback on screen
prevnumsubmitted = 0;
prevnumcollected = 0;
prevnumbusy      = 0;

% determine the initial job order, small numbers are submitted first
if strcmp(order, 'random')
  priority = randperm(numjob);
elseif strcmp(order, 'original')
  priority = 1:numjob;
else
  error('unsupported order');
end

% post all jobs and gather their results
while ~all(submitted) || ~all(collected)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 1: submit the jobs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % select all jobs that still need to be submitted
  submit = find(~submitted);

  if ~isempty(submit) && (sum(submitted)-sum(collected))<MaxBusy

    % determine the job to submit, the one with the smallest priority number goes first
    [dum, sel] = min(priority(submit));
    submit = submit(sel(1));

    % redistribute the input arguments
    argin = cell(1, numargin);
    for j=1:numargin
      argin{j} = varargin{j}{submit};
    end

    % submit the job for execution
    ws = warning('off', 'FieldTrip:peer:noSlaveAvailable');
    % peerfeval will give a warning if the submission timed out
    [curjobid curputtime] = peerfeval(fname, argin{:}, 'timeout', 5, 'memreq', memreq, 'timreq', timreq, 'diary', diary);
    warning(ws);

    if ~isempty(curjobid)
      % fprintf('submitted job %d\n', submit);
      jobid(submit)      = curjobid;
      puttime(submit)    = curputtime;
      submitted(submit)  = true;
      submittime(submit) = toc(stopwatch);
      clear curjobid curputtime
    end

    clear argin
  end % if ~isempty(submitlist)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 2: collect the job results that have finished sofar
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % get the list of job results
  joblist = peer('joblist');

  % get the output of all jobs that have finished
  for i=1:numel(joblist)

    % figure out to which job these results belong
    collect = find(jobid == joblist(i).jobid);

    if isempty(collect) && ~isempty(resubmitted)
      % it might be that these results are from a previously resubmitted job
      collect = [resubmitted([resubmitted.jobid] == joblist(i).jobid).jobnum];
      if ~isempty(collect) && ~collected(collect)
        % forget the resubmitted job, take these results instead
        warning('the original job %d did return, reverting to its original results', collect);
      end
    end

    if isempty(collect)
      % this job is not interesting to collect, probably it reflects some junk
      % from a previous call or a failed resubmission
      peer('clear', joblist(i).jobid);
      continue;
    end

    if collected(collect)
      % this job is the result of a resubmission, where the original result did return
      peer('clear', joblist(i).jobid);
      continue;
    end

    % collect the output arguments
    [argout, options] = peerget(joblist(i).jobid, 'timeout', inf, 'output', 'cell', 'diary', diary, 'StopOnError', StopOnError);

    % fprintf('collected job %d\n', collect);
    collected(collect)   = true;
    collecttime(collect) = toc(stopwatch);

    % redistribute the output arguments
    for j=1:numargout
      varargout{j}{collect} = argout{j};
    end

    % gather the job statistics
    % these are empty in case an error happened during remote evaluation, therefore the default value of NaN is specified
    timused(collect) = keyval('timused', options, nan);
    memused(collect) = keyval('memused', options, nan);

    prev_timreq = timreq;
    prev_memreq = memreq;

    if any(collected)
      % update the estimate of the time and memory that will be needed for the next job
      timreq = nanmax(timused);
      memreq = nanmax(memused);
    elseif ~any(collected) && any(submitted) && any(busy)
      % update based on the time spent waiting sofar for the first job to return
      elapsed = toc(stopwatch) - min(submittime(submitted(busy)));
      timreq  = max(timreq, elapsed);
    end

    % give some feedback
    if memreq~=prev_memreq
      memreq_in_mb = memreq/(1024*1024);
      if memreq_in_mb<100
        fprintf('updating memreq to %.3f MB\n', memreq_in_mb);
      else
        fprintf('updating memreq to %.0f MB\n', memreq_in_mb);
      end
    end

    % give some feedback
    if timreq~=prev_timreq
      if timreq<100
        fprintf('updating timreq to %.3f s\n', timreq);
      else
        fprintf('updating timreq to %.0f s\n', timreq);
      end
    end

  end % for joblist

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 3: flag jobs that take too long for resubmission
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % search for jobs that were submitted but that are still not busy after 30 seconds
  % this happens if the peerslave is not able to get a matlab license
  elapsed = toc(stopwatch) - submittime;
  elapsed(~submitted) = 0;
  elapsed(collected)  = 0;
  elapsed(busy)       = 0;
  sel = find(elapsed>30);

  for i=1:length(sel)
    warning('resubmitting job %d because it takes too long to get started', sel(i));
    % remember the job that will be resubmitted, it still might return its results
    resubmitted(end+1).jobnum = sel(i);
    resubmitted(end  ).jobid  = jobid(sel(i));
    resubmitted(end  ).time   = toc(stopwatch);
    resubmitted(end  ).reason = 'startup';

    % reset all job information, this will cause it to be automatically resubmitted
    jobid      (sel(i)) = nan;
    puttime    (sel(i)) = nan;
    timused    (sel(i)) = nan;
    memused    (sel(i)) = nan;
    submitted  (sel(i)) = false;
    collected  (sel(i)) = false;
    busy       (sel(i)) = false;
    submittime (sel(i)) = inf;
    collecttime(sel(i)) = inf;

    % increase the priority number, the resubmission should as late as possible
    % to increase the chance of the original job returning its results
    priority(sel(i)) = max(priority)+1;
  end

  % search for jobs that take too long to return their results
  % use an estimate of the time it requires a job to complete

  if ~isempty(ResubmitTime)
    % use the user-specified amount
    estimated = ResubmitTime;
  elseif any(collected)
    % estimate the time that it took the collected jobs to finish
    estimated_min = min(collecttime(collected) - submittime(collected));
    estimated_max = max(collecttime(collected) - submittime(collected));
    estimated_avg = max(collecttime(collected) - submittime(collected));
    % estimate the expected time of the jobs, assuming a "normal" distribution
    % the rationale for the estimate is the average plus N times the standard deviation
    estimated = estimated_avg + 2*(estimated_max - estimated_min);
  elseif ~isempty(timreq)
    % assume that it will not take more than 3x the required time
    estimated = 3*timreq;
  else
    % it is not possible to estimate the time that a job will take
    estimated = inf;
  end

  % add some time to allow the matlab engine to start
  estimated = estimated + 30;

  % test whether one of the submitted jobs should be resubmitted
  elapsed = toc(stopwatch) - submittime;
  sel = find(submitted & ~collected & (elapsed>estimated));

  for i=1:length(sel)
    warning('resubmitting job %d because it takes too long to finish (estimated = %f, elapsed = %f)', sel(i), estimated, elapsed(sel(i)));
    % remember the job that will be resubmitted, it still might return its results
    resubmitted(end+1).jobnum = sel(i);
    resubmitted(end  ).jobid  = jobid(sel(i));
    resubmitted(end  ).time   = toc(stopwatch);
    resubmitted(end  ).reason = 'duration';

    % reset all job information, this will cause it to be automatically resubmitted
    jobid      (sel(i)) = nan;
    puttime    (sel(i)) = nan;
    timused    (sel(i)) = nan;
    memused    (sel(i)) = nan;
    submitted  (sel(i)) = false;
    collected  (sel(i)) = false;
    busy       (sel(i)) = false;
    submittime (sel(i)) = inf;
    collecttime(sel(i)) = inf;

    % increase the priority number, the resubmission should as late as possible
    % to increase the chance of the original job returning its results
    priority(sel(i)) = max(priority)+1;
  end

  if ~all(collected)
    % wait a little bit, then try again to submit or collect a job
    pause(sleep);
  end

  % get the list of jobs that are busy
  busylist = peerlist('busy');
  busy(:)  = false;
  if ~isempty(busylist)
    current    = [busylist.current];
    [dum, sel] = intersect(jobid, [current.jobid]);
    busy(sel)  = true; % this indicates that the job execution is currently busy
  end

  if sum(collected)>prevnumcollected || sum(busy)~=prevnumbusy
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d, busy %d, speedup %.1f\n', sum(submitted), numel(submitted), sum(collected), numel(collected), sum(busy), nansum(timused(collected))/toc(stopwatch));
  end

  prevnumsubmitted = sum(submitted);
  prevnumcollected = sum(collected);
  prevnumbusy      = sum(busy);

end % while not all jobs have finished

if numargout>0 && UniformOutput
  % check whether the output can be converted to a uniform one
  for i=1:numel(varargout)
    for j=1:numel(varargout{i})
      if numel(varargout{i}{j})~=1
        % this error message is consistent with the one from cellfun
        error('Non-scalar in Uniform output, at index %d, output %d. Set ''UniformOutput'' to false.', j, i);
      end
    end
  end

  % convert the output to a uniform one
  for i=1:numargout
    varargout{i} = [varargout{i}{:}];
  end
end

% compare the time used inside this function with the total execution time
fprintf('computational time = %.1f sec, elapsed = %.1f sec, speedup %.1f x (memreq = %s, timreq = %s)\n', nansum(timused), toc(stopwatch), nansum(timused)/toc(stopwatch), print_mem(memreq), print_tim(timreq));

if all(puttime>timused)
  % FIXME this could be detected in the loop above, and the strategy could automatically
  % be adjusted from using the peers to local execution
  warning('copying the jobs over the network took more time than their execution');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this masks the regular version, this one only updates the status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function peerzombie
peer('status', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nanmax(x)
y = max(x(~isnan(x(:))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nanmin(x)
y = min(x(~isnan(x(:))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nanmean(x)
x = x(~isnan(x(:)));
y = mean(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nanstd(x)
x = x(~isnan(x(:)));
y = std(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nansum(x)
x = x(~isnan(x(:)));
y = sum(x);

