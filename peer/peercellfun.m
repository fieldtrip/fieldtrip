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
%   memreq         = number
%   timreq         = number
%   sleep          = number
%   diary          = string, can be 'always', 'never', 'warning', 'error' (default = 'error')
%   timcv          = coefficient of variation of the time required for the jobs (default is automatic)
%
% Example
%   fname = 'power';
%   x1    = {1, 2, 3, 4, 5};
%   x2    = {2, 2, 2, 2, 2};
%   y     = peercellfun(fname, x1, x2);
%
% See also CELLFUN, PEERMASTER, PEERFEVAL, PEERLIST, PEERINFO

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

% get the optional input arguments
UniformOutput = keyval('UniformOutput', optarg); if isempty(UniformOutput), UniformOutput = true; end
memreq  = keyval('memreq',  optarg); if isempty(memreq), memreq=1024^3;         end % assume 1 GB
timreq  = keyval('timreq',  optarg); if isempty(timreq), timreq=3600;           end % assume 1 hour
sleep   = keyval('sleep',   optarg); if isempty(sleep),  sleep=0.05;            end
diary   = keyval('diary',   optarg); if isempty(diary),  diary='error';         end
timcv   = keyval('timcv',   optarg); % default is empty, which will cause the range to be estimated

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

% post all jobs and gather their results
while ~all(submitted) || ~all(collected)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 1: submit the jobs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % select all jobs that still need to be submitted
  submit = find(~submitted);

  if ~isempty(submit)

    % pick the first job, or alternatively pick a random job
    % pick = 1;
    pick = ceil(rand(1)*length(submit));
    submit = submit(pick);

    if any(collected)
      prev_timreq = timreq;
      prev_memreq = memreq;
      % update the estimate of the time and memory that will be needed for the next job
      timreq = nanmax(timused);
      memreq = nanmax(memused);
      if timreq~=prev_timreq
        if timreq<100
          fprintf('updating timreq to %.3f s\n', timreq);
        else
          fprintf('updating timreq to %.0f s\n', timreq);
        end
      end
      if memreq~=prev_memreq
        memreq_in_mb = memreq/(1024*1024);
        if memreq_in_mb<100
          fprintf('updating memreq to %.3f MB\n', memreq_in_mb);
        else
          fprintf('updating memreq to %.0f MB\n', memreq_in_mb);
        end
      end
    elseif ~any(collected) && any(submitted)
      prev_timreq = timreq;
      % update based on the time spent waiting sofar for the first job to return
      elapsed = toc(stopwatch) - min(submittime(submitted));
      timreq  = max(timreq, elapsed);
      if timreq~=prev_timreq
        fprintf('updating timreq to %d\n', timreq);
      end
    end

    % redistribute the input arguments
    argin = cell(1, numargin);
    for j=1:numargin
      argin{j} = varargin{j}{submit};
    end

    % get a list of the busy slaves, used for feedback in case peerfeval times out
    busy = peerlist('busy');

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
    else
      if ~isempty(busy)
        % select only the slaves that are busy with your jobs
        current = [busy.current];
        info = peerinfo;
        busy = busy([current.hostid]==info.hostid);
        clear current info
      end
      if isempty(busy)
        warning('none of the slaves seems to be busy with your jobs');
      end
    end

    clear argin
  end % if ~isempty(submitlist)

  if sum(submitted)>prevnumsubmitted
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d, busy %d, speedup %.1f\n', sum(submitted), numel(submitted), sum(collected), numel(collected), sum(submitted)-sum(collected), sum(timused(collected))/toc(stopwatch));
  end

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
      if ~isempty(collect)
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

    % collect the output arguments
    [argout, options] = peerget(joblist(i).jobid, 'timeout', inf, 'output', 'cell', 'diary', diary);

    % fprintf('collected job %d\n', collect);
    collected(collect)   = true;
    collecttime(collect) = toc(stopwatch);

    % redistribute the output arguments
    for j=1:numargout
      varargout{j}{collect} = argout{j};
    end

    % gather the job statistics
    timused(collect) = keyval('timused', options);
    memused(collect) = keyval('memused', options);

  end % for joblist

  if sum(collected)>prevnumcollected
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d, busy %d, speedup %.1f\n', sum(submitted), numel(submitted), sum(collected), numel(collected), sum(submitted)-sum(collected), sum(timused(collected))/toc(stopwatch));
  end

  prevnumsubmitted = sum(submitted);
  prevnumcollected = sum(collected);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 3: flag jobs that take too long for resubmission
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % check for jobs that are taking too long to finish
  if all(submitted) && any(collected) && ~all(collected)

    % estimate the elapsed time for all jobs
    elapsed = toc(stopwatch) - submittime;

    % estimate the time that it took the collected jobs to finish
    estimated_min = min(collecttime(collected) - submittime(collected));
    estimated_max = max(collecttime(collected) - submittime(collected));
    estimated_avg = estimated_max; % the maximum is used instead of the mean

    % the rationale for the estimate is the mean plus 2x the standard deviation
    if isempty(timcv)
      % instead of the standard deviation the min-max range (divided by two) is used
      estimated = estimated_avg + (estimated_max - estimated_min);
      % take into account that the estimate is inaccurate in case of few collected jobs
      estimated = estimated * (1 + 1/(1+log10(sum(collected))));
    else
      % the coefficient of variation (CV) is a normalized measure of dispersion of a distribution
      % it is defined as the ratio of the standard deviation to the mean
      estimated = estimated_avg + 2*timcv*estimated_avg;
    end

    % test whether one of the submitted jobs should be resubmitted
    sel = find(submitted & ~collected & (elapsed>estimated));

    for i=1:length(sel)
      warning('resubmitting job %d because it takes too long to finish (estimated = %f, elapsed = %f)', sel(i), estimated, elapsed(sel(i)));
      % remember the job that will be resubmitted, it still might return its results
      resubmitted(end+1).jobnum = sel(i);
      resubmitted(end  ).jobid  = jobid(sel(i));
      % reset all job information, this will cause it to be automatically resubmitted
      jobid      (sel(i)) = nan;
      puttime    (sel(i)) = nan;
      timused    (sel(i)) = nan;
      memused    (sel(i)) = nan;
      submitted  (sel(i)) = false;
      collected  (sel(i)) = false;
      submittime (sel(i)) = inf;
      collecttime(sel(i)) = inf;
    end
  end % resubmitting

  if ~all(collected)
    % wait a little bit, then try again to submit or collect a job
    pause(sleep);
  end

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
fprintf('computational time = %.1f sec, elapsed time = %.1f sec, approximate speedup %.1f x\n', sum(timused), toc(stopwatch), sum(timused)/toc(stopwatch));

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

