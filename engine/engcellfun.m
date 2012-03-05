function varargout = engcellfun(fname, varargin)

% ENGCELLFUN applies a function to each element of a cell-array. The function
% execution is done in parallel on locally or remotely running MATLAB engines.
%
% Use as
%   argout = engcellfun(fname, x1, x2, ...)
%
% Optional arguments can be specified in key-value pairs and can include
%   UniformOutput  = boolean (default = false)
%   StopOnError    = boolean (default = true)
%   diary          = string, can be 'always', 'never', 'warning', 'error' (default = 'error')
%   order          = string, can be 'random' or 'original' (default = 'random')
%
%  Example
%    x1 = {1, 2, 3, 4, 5};
%    x2 = {2, 2, 2, 2, 2};
%    engpool open 4
%    y  = engcellfun(@power, x1, x2);
%    engpool close
%
% See also ENGPOOL, ENGFEVAL, ENGGET

% -----------------------------------------------------------------------
% Copyright (C) 2012, Robert Oostenveld
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

if matlabversion(7.8, Inf)
  % switch to zombie when finished or when ctrl-c gets pressed
  % the onCleanup function does not exist for older versions
  onCleanup(@cleanupfun);
end

% locate the begin of the optional key-value arguments
optbeg = find(cellfun(@ischar, varargin));
optarg = varargin(optbeg:end);

% get the optional input arguments
UniformOutput = ft_getopt(optarg, 'UniformOutput', false   );
StopOnError   = ft_getopt(optarg, 'StopOnError',   true    );
diary         = ft_getopt(optarg, 'diary',         'error' );   % 'always', 'never', 'warning', 'error'
order         = ft_getopt(optarg, 'order',         'original'); % 'random', 'original'

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

% check the availability of the engines
pool = engpool('info');
pool   = engpool('info');
isbusy = false(1,numel(pool));
hasjob = false(1,numel(pool));
for i=1:numel(pool)
  isbusy(i) = engine('isbusy', i);
  hasjob(i) = ~isempty(pool{i});
end

if isempty(pool)
  error('the engine pool is empty');
end

if all(isbusy)
  warning('there is no engine available, reverting to local cellfun');
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

% remove any remains from an aborted previous call
cleanupfun

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
  
  pool   = engpool('info');
  isbusy = false(1,numel(pool));
  hasjob = false(1,numel(pool));
  for i=1:numel(pool)
    isbusy(i) = engine('isbusy', i);
    hasjob(i) = ~isempty(pool{i});
  end
  
  if ~isempty(submit) && ~all(hasjob | isbusy)
    
    % determine the job to submit, the one with the smallest priority number goes first
    [dum, sel] = min(priority(submit));
    submit = submit(sel(1));
    
    % redistribute the input arguments
    argin = cell(1, numargin);
    for j=1:numargin
      argin{j} = varargin{j}{submit};
    end
    
    % submit the job for execution
    [curjobid curputtime] = engfeval(fname, argin{:}, 'diary', diary);
    
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
  
  % check the availability of the engines
  pool   = engpool('info');
  isbusy = false(1,numel(pool));
  hasjob = false(1,numel(pool));
  for i=1:numel(pool)
    isbusy(i) = engine('isbusy', i);
    hasjob(i) = ~isempty(pool{i});
  end
  
  if sum(collected)>prevnumcollected || sum(isbusy)~=prevnumbusy
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d, busy %d, speedup %.1f\n', sum(submitted), numel(submitted), sum(collected), numel(collected), sum(isbusy), nansum(timused(collected))/toc(stopwatch));
  end
  
  prevnumcollected = sum(collected);
  prevnumbusy      = sum(isbusy);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 2: collect the job results that have finished sofar
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % check whether one of the engines finished
  ready = hasjob & ~isbusy;
  
  % get the output of that engine
  if any(ready)
    
    % select the first engine that is ready
    ready = find(ready, 1, 'first');
    % figure out to which job its results belong
    collect = find(jobid == pool{ready});
    
    if isempty(collect)
      % this job is not interesting to collect, probably it reflects some junk
      % from a previous call or a failed resubmission
      keyboard
      % FIXME
      continue;
    end
    
    if collected(collect)
      % this job is the result of a resubmission, where the original result did return
      keyboard
      % FIXME
      continue;
    end
    
    % collect the output arguments
    try
      ws = warning('Off','Backtrace');
      [argout, options] = engget(pool{ready}, 'output', 'cell', 'diary', diary, 'StopOnError', StopOnError);
      warning(ws);
    catch
      keyboard
      % FIXME
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
    
  end % collect the selected job
  
  if sum(collected)>prevnumcollected
    % give an update of the progress
    fprintf('submitted %d/%d, collected %d/%d, busy %d, speedup %.1f\n', sum(submitted), numel(submitted), sum(collected), numel(collected), sum(isbusy), nansum(timused(collected))/toc(stopwatch));
  end
  
  if all(isbusy)
    % wait a little bit for the running jobs to complete
    pause(0.1);
  end
  
  prevnumcollected = sum(collected);
  prevnumbusy      = sum(isbusy);
  
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
fprintf('computational time = %.1f sec, elapsed = %.1f sec, speedup %.1f x\n', nansum(timused), toc(stopwatch), nansum(timused)/toc(stopwatch));

if all(puttime>timused)
  % FIXME this could be detected in the loop above, and the strategy could automatically
  % be adjusted from using the engines to local execution
  warning('copying the jobs over the network took more time than their execution');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanupfun
pool = engpool('info');
for i=1:length(pool)
  engpool('release', i);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nansum(x)
x = x(~isnan(x(:)));
y = sum(x);
