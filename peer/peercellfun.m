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
% See also CELLFUN, PEERMASTER, PEERFEVAL, PEERGET

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

  % select one of the jobs to be submitted
  submit = find(~submitted);                % select all jobs that still need to be submitted
  submit = submit(randperm(numel(submit))); % randomize the order of the list of jobs

  if ~isempty(submit)
    % take the first job from the list
    submit = submit(1);

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

    % submit the job for execution
    [jobid(submit) puttime(submit)] = peerfeval(fname, argin{:}, 'timeout', inf, 'memreq', memreq, 'timreq', timreq, 'diary', diary);
    % fprintf('submitted job %d\n', submit);
    submitted(submit)  = true;
    submittime(submit) = toc(stopwatch);

    clear argin
  end % if ~isempty(submit)

  if sum(submitted)>prevnumsubmitted
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d, busy %d\n', sum(submitted), numel(submitted), sum(collected), numel(collected), sum(submitted)-sum(collected));
  end

  joblist = peer('joblist');

  % get the results of all jobs that have finished
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
    fprintf('submitted %d/%d, collected %d/%d, busy %d\n', sum(submitted), numel(submitted), sum(collected), numel(collected), sum(submitted)-sum(collected));
  end

  prevnumsubmitted = sum(submitted);
  prevnumcollected = sum(collected);

  % check for jobs that are taking too long to finish
  if all(submitted) && any(collected) && ~all(collected)
    % test whether one of the jobs should be resubmitted
    sel = find(~collected, 1);
    elapsed = toc(stopwatch) - submittime(sel);

    % estimate the time that it took the other jobs to finish
    estimated_min = min(collecttime(collected) - submittime(collected));
    estimated_max = max(collecttime(collected) - submittime(collected));
    estimated_avg = estimated_max; % the maximum is used instead of the mean

    % the rationale for the estimate is the mean plus 2x the standard deviation
    if isempty(timcv)
      % instead of the standard deviation the min-max range (divided by two) is used
      estimated = estimated_avg + (estimated_max - estimated_min);
      % take into account that the estimate is inaccurate in case of few collected jobs
      estimated = (1 + 2^(-sum(collected))) * estimated;
    else
      % the coefficient of variation (CV) is a normalized measure of dispersion of a distribution
      % it is defined as the ratio of the standard deviation to the mean
      estimated = estimated_avg + 2*timcv*estimated_avg;
    end

    if elapsed>estimated
      warning('resubmitting job %d because it takes too long to finish (estimated = %f, elapsed = %f)', sel, estimated, elapsed);
      % remember the job that will be resubmitted, it still might return its results
      resubmitted(end+1).jobnum = sel;
      resubmitted(end  ).jobid  = jobid(sel);
      % reset all job information, this will cause it to be automatically resubmitted
      jobid      (sel) = nan;
      puttime    (sel) = nan;
      timused    (sel) = nan;
      memused    (sel) = nan;
      submitted  (sel) = false;
      collected  (sel) = false;
      submittime (sel) = inf;
      collecttime(sel) = inf;
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
fprintf('approximate speedup ratio %f\n', sum(timused)/toc(stopwatch));

if all(puttime>timused)
  % FIXME this could be detected in the loop above, and the strategy could automatically
  % be adjusted from using the peers to local execution
  warning('copying the jobs over the network took more time than their execution');
end

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

