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
%   diary          = string, can be 'always', 'warning', 'error' (default = 'error')
%
% Example
%   fname = 'power';
%   x1    = {1, 2, 3, 4, 5};
%   x2    = {2, 2, 2, 2, 2};
%   y     = peercellfun(fname, x1, x2);
%
% See also CELLFUN, PEERMASTER, PEERFEVAL, PEERGET
%
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
memreq  = keyval('memreq',  optarg); if isempty(memreq), memreq=0;       end
timreq  = keyval('timreq',  optarg); if isempty(timreq), timreq=3600;    end % assume that it will take one hour
sleep   = keyval('sleep',   optarg); if isempty(sleep),  sleep=0.01;     end
diary   = keyval('diary',   optarg); if isempty(diary),  diary='error';  end

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
catch
  me = lasterror;
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

% start the timer
stopwatch = tic;

% these are used for printing feedback on screen
prevnumsubmitted = 0;
prevnumcollected = 0;

% post all jobs and gather their results
while ~all(submitted) || ~all(collected)

  % select one of the jobs to be submitted
  submit = find(~submitted);                % select all jobs that still need to be submitted
  submit = submit(randperm(numel(submit))); % change into a random order

  if ~isempty(submit)
    % take the first job from the random list
    submit = submit(1);

    if any(collected)
      % update the estimate of the time and memory that will be needed for the next job
      timreq = nanmax(timused);
      memreq = nanmax(memused);
    elseif ~any(collected) && any(submitted)
      % update based on the time spent waiting sofar for the first job to return
      elapsed = toc(stopwatch) - min(submittime(submitted));
      timreq  = max(timreq, elapsed);
    end

    % redistribute the input arguments
    argin = cell(1, numargin);
    for j=1:numargin
      argin{j} = varargin{j}{submit};
    end

    [jobid(submit) puttime(submit)] = peerfeval(fname, argin{:}, 'timeout', inf, 'memreq', memreq, 'timreq', timreq);
    % fprintf('submitted job %d\n', submit);
    submitted(submit)  = true;
    submittime(submit) = toc(stopwatch);
    clear argin
  end

  if sum(submitted)>prevnumsubmitted
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d\n', sum(submitted), numel(submitted), sum(collected), numel(collected));
  end

  joblist = peer('joblist');

  % get the results of all jobs that have finished
  for i=1:numel(joblist)
    collect = find(jobid == joblist(i).jobid);
    if isempty(collect)
      % there must be some junk in the buffer from an aborted previous call
      [argout, options] = peerget(joblist(i).jobid, 'timeout', inf, 'output', 'cell', 'diary', diary);
    else
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
    end
  end

  if sum(collected)>prevnumcollected
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d\n', sum(submitted), numel(submitted), sum(collected), numel(collected));
  end

  if all(submitted) && ~all(collected)
    % wait a little bit and try to collect another job
    pause(sleep);
  end

  prevnumsubmitted = sum(submitted);
  prevnumcollected = sum(collected);

  % check for jobs that are taking too long to finish
  if all(submitted) && ~all(collected)
    % test whether one of the jobs should be resubmitted
    sel = find(~collected, 1);
    elapsed = toc(stopwatch) - submittime(sel);
    % estimate the time that it took the other jobs to finish
    estimated_min = min(collecttime(collected) - submittime(collected));
    estimated_max = max(collecttime(collected) - submittime(collected));
    % the rationale for the estimate is the mean plus 2x the standard deviation
    % except that instead of the mean the maximum is used and instead of the 
    % standard deviation the min-max range is used
    estimated     = estimated_max + (estimated_max - estimated_min);
    if elapsed>estimated
      warning('resubmitting job %d because it took too long to finish (estimated = %f, elapsed = %f)', sel, estimated, elapsed);
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

