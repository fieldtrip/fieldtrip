function varargout = qsubcellfun(fname, varargin)

% QSUBCELLFUN applies a function to each element of a cell-array. The
% function execution is done in parallel using the Torque, SGE or PBS 
% batch queue system.
%
% Use as
%   argout = qsubcellfun(fname, x1, x2, ...)
%
% This function has a number of optional arguments that have to passed
% as key-value pairs at the end of the list of input arguments. All other
% input arguments (including other key-value pairs) will be passed to the
% function to be evaluated.
%   UniformOutput  = boolean (default = false)
%   StopOnError    = boolean (default = true)
%   diary          = string, can be 'always', 'never', 'warning', 'error' (default = 'error')
%   timreq         = number, the time in seconds required to run a single job
%   memreq         = number, the memory in bytes required to run a single job
%
% It is required to give an estimate of the time and memory requirements of
% the individual jobs. The memory requirement of the MATLAB executible
% itself will automatically be added, just as the time requirement to start
% up a new MATLAB process. If you don't know what the memory and time
% requirements of your job are, you can get an estimate for them using
% TIC/TOC and MEMTIC/MEMTOC around a single execution of one of the jobs in
% your interactive MATLAB session. You can also start with very large
% estimates, e.g. 4*1024^3 bytes for the memory (which is 4GB) and 28800
% seconds for the time (which is 8 hours) and then run a single job through
% qsubcellfun. When the job returns, it will print the memory and time it
% required.
%
% Example
%   fname = 'power';
%   x1    = {1, 2, 3, 4, 5};
%   x2    = {2, 2, 2, 2, 2};
%   y     = qsubcellfun(fname, x1, x2, 'memreq', 1024^3, 'timreq', 60);
%
% In case you abort your call to qsubcellfun by pressing ctrl-C, the
% already submitted jobs will continue to run. You will also notice that
% a lot of temporary files remain in your working directory. To check the
% status and healthy execution of the jobs in the batch queuing system,
% you can use
%   qstat
%   qstat -an1
%   qstat -Q
% comands on the linux command line. To delete jobs from the queue  and
% to abort already running jobs, you can use
%   qdel <jobnumber>
%   qdel all
%
% See also QSUBFEVAL, CELLFUN, PEERCELLFUN, FEVAL, DFEVAL, DFEVALASYNC

% -----------------------------------------------------------------------
% Copyright (C) 2011, Robert Oostenveld
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

% this persistent variable is used to assign unique names to the jobs in each subsequent batch
persistent batch;

if matlabversion(7.8, Inf)
  % switch to zombie when finished or when ctrl-c gets pressed
  % the onCleanup function does not exist for older versions
  onCleanup(@cleanupfun);
end

% remove the persistent lists with job and pbs identifiers
clear qsublist

stopwatch = tic;

% locate the begin of the optional key-value arguments
optbeg = find(cellfun(@ischar, varargin));
optarg = varargin(optbeg:end);

% get the optional input arguments
UniformOutput = ft_getopt(optarg, 'UniformOutput', false   );
StopOnError   = ft_getopt(optarg, 'StopOnError',   true    );
diary         = ft_getopt(optarg, 'diary',         'error' ); % 'always', 'never', 'warning', 'error'
timreq        = ft_getopt(optarg, 'timreq');
memreq        = ft_getopt(optarg, 'memreq'); 

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

% prepare some arrays that are used for bookkeeping
jobid       = cell(1, numjob);
puttime     = nan(1, numjob);
timused     = nan(1, numjob);
memused     = nan(1, numjob);
submitted   = false(1, numjob);
collected   = false(1, numjob);
submittime  = inf(1, numjob);
collecttime = inf(1, numjob);

% keep a record of how many calls to this function have been made (useful for creating meaningful job IDs)
% subsequent batches of jobs will have different names
if isempty(batch)
  batch = 1;
else
  batch = batch+1;
end

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
    error('input argument #%d shoudl be a cell-array', i+1);
  end
  if numel(varargin{i})~=numjob
    error('inconsistent number of elements in input #%d', i+1);
  end
end

for submit=1:numjob
  % redistribute the input arguments
  argin = cell(1, numargin);
  for j=1:numargin
    argin{j} = varargin{j}{submit};
  end

  % submit the job
  [curjobid curputtime] = qsubfeval(fname, argin{:}, 'memreq', memreq, 'timreq', timreq, 'diary', diary, 'batch', batch);

  % fprintf('submitted job %d\n', submit);
  jobid{submit}      = curjobid;
  puttime(submit)    = curputtime;
  submitted(submit)  = true;
  submittime(submit) = toc(stopwatch);
  clear curjobid curputtime
end % for

while (~all(collected))
  % try to collect the jobs that have finished

  for collect=find(~collected)
    % this will return empty arguments if the job has not finished
    ws = warning('off', 'FieldTrip:qsub:jobNotAvailable');
    [argout, options] = qsubget(jobid{collect}, 'output', 'cell', 'diary', diary, 'StopOnError', StopOnError);
    warning(ws);

    if ~isempty(argout) || ~isempty(options)
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

    end  % if
  end % for

  pause(0.1);
end % while

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nansum(x)
x = x(~isnan(x(:)));
y = sum(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanupfun
% the qsublist function maintains a persistent list with all jobs
% request it to kill all the jobs and to cleanup all the files
qsublist('killall');

