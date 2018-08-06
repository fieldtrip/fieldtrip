function varargout = peercellfun(fname, varargin)

% PEERCELLFUN applies a function to each element of a cell-array. The
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
%   RetryOnError   = number, number of retries for failed jobs expressed as ratio (default = 0.05)
%   MaxBusy        = number, amount of slaves allowed to be busy (default = inf)
%   diary          = string, can be 'always', 'never', 'warning', 'error' (default = 'error')
%   timreq         = number, initial estimate for the time required to run a single job (default = 3600)
%   mintimreq      = number, minimum time required to run a single job (default is automatic)
%   memreq         = number, initial estimate for the memory required to run a single job (default = 2*1024^3)
%   minmemreq      = number, minimum memory required to run a single job (default is automatic)
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
%
% $Id$
% -----------------------------------------------------------------------

if ft_platform_supports('onCleanup')
  % switch to zombie when finished or when ctrl-c gets pressed
  % the onCleanup function does not exist for older versions
  onCleanup(@peerzombie);
end

% locate the begin of the optional key-value arguments
optbeg = find(cellfun(@ischar, varargin));
optarg = varargin(optbeg:end);

% get the optional input arguments
UniformOutput = ft_getopt(optarg, 'UniformOutput', false   );
StopOnError   = ft_getopt(optarg, 'StopOnError',   true    );
MaxBusy       = ft_getopt(optarg, 'MaxBusy',       inf     );
RetryOnError  = ft_getopt(optarg, 'RetryOnError',  0.050   ); % ratio, fraction of the total jobs
sleep         = ft_getopt(optarg, 'sleep',         0.050   ); % time in seconds
diary         = ft_getopt(optarg, 'diary',         'error' ); % 'always', 'never', 'warning', 'error'
order         = ft_getopt(optarg, 'order',         'random'); % 'random', 'original'
timreq        = ft_getopt(optarg, 'timreq',        []      ); 
mintimreq     = ft_getopt(optarg, 'mintimreq',     []      ); 
memreq        = ft_getopt(optarg, 'memreq',        []      ); % see below
minmemreq     = ft_getopt(optarg, 'minmemreq',     []      ); % see below

if isempty(timreq) && isempty(mintimreq)
  % assume an initial job duration of 1 hour
  % the time required by the jobs will be estimated and timreq will be auto-adjusted
  timreq    = 3600;
  mintimreq = 0;
elseif isempty(timreq)
  % use the specified mimimum as the initial value that a job required
  % it will be auto-adjusted to larger values, not to smaller values
  timreq    = mintimreq;
elseif isempty(mintimreq)
  % jobs will be killed by the slave if they take more than 3 times the estimated time at submission
  % use the user-supplied initial value, the minimum should not be less than 1/3 of that
  mintimreq = timreq/3;
end

if isempty(memreq) && isempty(minmemreq)
  % assume an initial memory requirement of 1 GB
  % the memory required by the jobs will be estimated and memreq will be auto-adjusted
  memreq    = 2*1024^3;
  minmemreq = 0;
elseif isempty(memreq)
  % use the specified mimimum as the initial value that a job required
  % it will be auto-adjusted to larger values, not to smaller values
  memreq    = minmemreq;
elseif isempty(minmemreq)
  % jobs will be killed by the slave if they take more than 1.5 times the estimated time at submission
  % use the user-supplied initial value, the minimum should not be less than 1/1.5 times that
  minmemreq = memreq/1.5;
end

% convert from 'yes'/'no' into boolean value
UniformOutput = istrue(UniformOutput);

% convert from a fraction into an integer number
RetryOnError = floor(RetryOnError * numel(varargin{1}));

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
catch
  % the "catch me" syntax is broken on MATLAB74, this fixes it
  nargout_err = lasterror;
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
elseif numargout>nargout
  % the number of output arguments is constrained by the users' call to this function
  numargout = nargout;
elseif nargout>numargout
  error('Too many output arguments.');
end

% check the input arguments
for i=1:numargin
  if ~isa(varargin{i}, 'cell')
    error('input argument #%d should be a cell-array', i+1);
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
lastseen    = inf(1, numjob);
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
prevtimreq       = timreq;
prevmemreq       = memreq;
  
  if any(collected)
    % update the estimate of the time and memory that will be needed for the next job
    % note that it cannot be updated if all collected jobs have failed (in case of stoponerror=false)
    if ~isempty(nanmax(timused))
      timreq = nanmax(timused);
      timreq = max(timreq, mintimreq);
      memreq = nanmax(memused);
      memreq = max(memreq, minmemreq);
    end
  end
  if any(submitted) && any(busy)
    % update based on the time already spent on the slowest job
    elapsed = toc(stopwatch) - min(submittime(submitted & busy));
    timreq  = max(timreq, elapsed);
    timreq  = max(timreq, mintimreq);
  end

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
    [curjobid curputtime] = peerfeval(fname, argin{:}, 'timeout', 5, 'memreq', memreq, 'timreq', timreq, 'diary', diary, 'nargout', numargout);
    warning(ws);

    if ~isempty(curjobid)
      % fprintf('submitted job %d\n', submit);
      jobid(submit)      = curjobid;
      puttime(submit)    = curputtime;
      submitted(submit)  = true;
      submittime(submit) = toc(stopwatch);
      clear curjobid curputtime

      % give some feedback
      if abs(memreq-prevmemreq)>1000
        fprintf('updating memreq to %s\n', print_mem(memreq));
      end

      % give some feedback
      if abs(timreq-prevtimreq)>1
        fprintf('updating timreq to %s\n', print_tim(timreq));
      end
    end

    clear argin
  end % if ~isempty(submitlist)

  % get the list of jobs that are busy
  busylist = peerlist('busy');
  busy(:)  = false;
  if ~isempty(busylist)
    current       = [busylist.current];
    [dum, sel]    = intersect(jobid, [current.jobid]);
    busy(sel)     = true;           % this indicates that the job execution is currently busy
    lastseen(sel) = toc(stopwatch); % keep track of when the job was seen the last time
  end

  if sum(collected)>prevnumcollected || sum(busy)~=prevnumbusy
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d, busy %d, speedup %.1f\n', sum(submitted), numel(submitted), sum(collected), numel(collected), sum(busy), nansum(timused(collected))/toc(stopwatch));
  end

  prevnumsubmitted = sum(submitted);
  prevnumcollected = sum(collected);
  prevnumbusy      = sum(busy);

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
    try
      ws = warning('Off','Backtrace');
      [argout, options] = peerget(joblist(i).jobid, 'timeout', inf, 'output', 'cell', 'diary', diary, 'StopOnError', StopOnError);
      warning(ws);
    catch
      % the "catch me" syntax is broken on MATLAB74, this fixes it
      peerget_err = lasterror;

      % the peerslave command line executable itself can return a number of errors
      %  1) could not start the MATLAB engine
      %  2) failed to execute the job (argin)
      %  3) failed to execute the job (optin)
      %  4) failed to execute the job (eval)
      %  5) failed to execute the job (argout)
      %  6) failed to execute the job (optout)
      %  7) failed to execute the job
      % errors 1-3 are not the users fault and happen prior to execution, therefore they should always result in a resubmission

      if ~isempty(strfind(peerget_err.message, 'could not start the MATLAB engine')) || ...
         ~isempty(strfind(peerget_err.message, 'failed to execute the job (argin)')) || ...
         ~isempty(strfind(peerget_err.message, 'failed to execute the job (optin)'))
        % this is due to a license problem or a memory problem
        if ~isempty(strfind(peerget_err.message, 'could not start the MATLAB engine'))
          warning('resubmitting job %d because the MATLAB engine could not get a license', collect);
        end
        % reset all job information, this will cause it to be automatically resubmitted
        jobid      (collect) = nan;
        puttime    (collect) = nan;
        timused    (collect) = nan;
        memused    (collect) = nan;
        submitted  (collect) = false;
        collected  (collect) = false;
        busy       (collect) = false;
        lastseen   (collect) = inf;
        submittime (collect) = inf;
        collecttime(collect) = inf;
        continue
      else
        % the returned error is more serious and requires the users attention
        if RetryOnError>0
          % dercease the counter for the remaining retries
          RetryOnError = RetryOnError - 1;
          % give the user some information
	      fprintf('an error was detected during the execution of job %d\n', collect);
          fprintf('??? %s\n', peerget_err.message);
	      fprintf('resubmitting the failed job (%d retries remaining)\n', RetryOnError);
          % reset all job information, this will cause it to be automatically resubmitted
          jobid      (collect) = nan;
          puttime    (collect) = nan;
          timused    (collect) = nan;
          memused    (collect) = nan;
          submitted  (collect) = false;
          collected  (collect) = false;
          busy       (collect) = false;
          lastseen   (collect) = inf;
          submittime (collect) = inf;
          collecttime(collect) = inf;
          continue
        else
	      fprintf('an error was detected during the execution of job %d\n', collect);
          rethrow(peerget_err);
        end
      end
    end
    
    % fprintf('collected job %d\n', collect);
    collected(collect)   = true;
    collecttime(collect) = toc(stopwatch);

    % redistribute the output arguments
    for j=1:numargout
      varargout{j}{collect} = argout{j};
    end

    % gather the job statistics
    % these are empty in case an error happened during remote evaluation, therefore the default value of NaN is specified
    timused(collect) = ft_getopt(options, 'timused', nan);
    memused(collect) = ft_getopt(options, 'memused', nan);

  end % for joblist

  prevtimreq = timreq;
  prevmemreq = memreq;
  
  if any(collected)
    % update the estimate of the time and memory that will be needed for the next job
    % note that it cannot be updated if all collected jobs have failed (in case of stoponerror=false)
    if ~isempty(nanmax(timused))
      timreq = nanmax(timused);
      timreq = max(timreq, mintimreq);
      memreq = nanmax(memused);
      memreq = max(memreq, minmemreq);
    end
  end
  if any(submitted) && any(busy)
    % update based on the time already spent on the slowest job
    elapsed = toc(stopwatch) - min(submittime(submitted & busy));
    timreq  = max(timreq, elapsed);
    timreq  = max(timreq, mintimreq);
  end

  % get the list of jobs that are busy
  busylist = peerlist('busy');
  busy(:)  = false;
  if ~isempty(busylist)
    current       = [busylist.current];
    [dum, sel]    = intersect(jobid, [current.jobid]);
    busy(sel)     = true;           % this indicates that the job execution is currently busy
    lastseen(sel) = toc(stopwatch); % keep track of when the job was seen the last time
  end

  if sum(collected)>prevnumcollected || sum(busy)~=prevnumbusy
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d, busy %d, speedup %.1f\n', sum(submitted), numel(submitted), sum(collected), numel(collected), sum(busy), nansum(timused(collected))/toc(stopwatch));
  end

  prevnumsubmitted = sum(submitted);
  prevnumcollected = sum(collected);
  prevnumbusy      = sum(busy);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 3: flag jobs that take too long for resubmission
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % this is only a warning, no action is taken here
  % sel = find((toc(stopwatch)-lastseen)>60);
  % for i=1:length(sel)
  %   warning('job %d has not been seen for 60 seconds\n', sel(i));
  % end

  % search for jobs that were submitted but that are still not busy after 60 seconds
  % this happens if the peerslave is not able to get a MATLAB license
  elapsed = toc(stopwatch) - submittime;
  elapsed(~submitted)       = 0;
  elapsed(collected)        = 0;
  elapsed(~isinf(lastseen)) = 0; % once started there is no reason to resubmit "because it takes too long to get started"
  sel = find(elapsed>60);

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
    lastseen   (sel(i)) = inf;
    submittime (sel(i)) = inf;
    collecttime(sel(i)) = inf;

    % increase the priority number, the resubmission should as late as possible
    % to increase the chance of the original job returning its results
    priority(sel(i)) = max(priority)+1;
  end

  % search for jobs that take too long to return their results
  % use an estimate of the time it requires a job to complete

  % assume that it will not take more than 3x the required time
  % this is also what is used by the peerslave to kill the job
  estimated = 3*timreq;

  % add some time to allow the MATLAB engine to start
  estimated = estimated + 60;

  % test whether one of the submitted jobs should be resubmitted
  elapsed = toc(stopwatch) - submittime;
  sel = find(submitted & ~collected & (elapsed>estimated));

  for i=1:length(sel)
    warning('resubmitting job %d because it takes too long to finish (estimated = %s)', sel(i), print_tim(estimated));
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
    lastseen   (sel(i)) = inf;
    submittime (sel(i)) = inf;
    collecttime(sel(i)) = inf;

    % increase the priority number, the resubmission should as late as possible
    % to increase the chance of the original job returning its results
    priority(sel(i)) = max(priority)+1;
  end

  if all(submitted)
    % wait a little bit, then try again to submit or collect a job
    pause(sleep);
  end

  % get the list of jobs that are busy
  busylist = peerlist('busy');
  busy(:)  = false;
  if ~isempty(busylist)
    current       = [busylist.current];
    [dum, sel]    = intersect(jobid, [current.jobid]);
    busy(sel)     = true;           % this indicates that the job execution is currently busy
    lastseen(sel) = toc(stopwatch); % keep track of when the job was seen the last time
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

% ensure the output to have the same size/dimensions as the input
for i=1:numargout
  varargout{i} = reshape(varargout{i}, size(varargin{1}));
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

