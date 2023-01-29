function varargout = qsubcellfun(fname, varargin)

% QSUBCELLFUN applies a function to each element of a cell-array. The
% function execution is done in parallel using the Torque, SGE, PBS or
% SLURM batch queue system.
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
%   timoverhead    = number in seconds, how much time to allow MATLAB to start (default = 180 seconds)
%   memreq         = number, the memory in bytes required to run a single job
%   memoverhead    = number in bytes, how much memory to account for MATLAB itself (default depends on the MATLAB version)
%   stack          = number, stack multiple jobs in a single qsub job (default = 'auto')
%   backend        = string, can be 'torque', 'sge', 'slurm', 'lsf', 'system', 'local' (default is automatic)
%   batchid        = string, to identify the jobs in the queue (default is user_host_pid_batch)
%   compile        = string, can be 'auto', 'yes', 'no' (default = 'no')
%   queue          = string, which queue to submit the job in (default is empty)
%   options        = string, additional options that will be passed to qsub/srun (default is empty)
%   matlabcmd      = string, the Linux command line to start MATLAB on the compute nodes (default is automatic
%   display        = 'yes' or 'no', whether the nodisplay option should be passed to MATLAB (default = 'no', meaning nodisplay)
%   jvm            = 'yes' or 'no', whether the nojvm option should be passed to MATLAB (default = 'yes', meaning with jvm)
%   rerunable      = 'yes' or 'no', whether the job can be restarted on a torque/maui/moab cluster (default = 'no')
%   sleep          = number, time in seconds to wait between checks for job completion (default = 0.5 s)
%
% It is required to give an estimate of the time and memory requirements of
% the individual jobs. The memory requirement of the MATLAB executable
% itself will automatically be added, just as the time required to start
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
%   y     = qsubcellfun(fname, x1, x2, 'memreq', 1024^3, 'timreq', 300);
%
% Using the compile=yes or compile=auto option, you can compile your
% function into a stand-alone executable that can be executed on the cluster
% without requiring additional MATLAB licenses. You can also call the
% QSUBCOMPILE function prior to calling QSUBCELLFUN. If you plan multiple
% batches of the same function, compiling it prior to QSUBCELLFUN is more
% efficient. In that case you will have to delete the compiled executable
% yourself once you are done.
%
% In case you abort your call to qsubcellfun by pressing ctrl-c,
% the already submitted jobs will be canceled. Some small temporary
% files might remain in your working directory.
%
% To check the the status and healthy execution of the jobs on the Torque
% batch queuing system, you can use
%   qstat
%   qstat -an1
%   qstat -Q
% comands on the linux command line. To delete jobs from the Torque batch
% queue and to abort already running jobs, you can use
%   qdel <jobnumber>
%   qdel all
%
% See also QSUBCOMPILE, QSUBFEVAL, CELLFUN, PEERCELLFUN, FEVAL, BATCH

% -----------------------------------------------------------------------
% Copyright (C) 2011-2016, Robert Oostenveld
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
  % switch to zombie when finished or when Ctrl-C gets pressed
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
timoverhead   = ft_getopt(optarg, 'timoverhead', []); % default is set in QSUBFEVAL
memoverhead   = ft_getopt(optarg, 'memoverhead', []); % default is set in QSUBFEVAL
stack         = ft_getopt(optarg, 'stack',   'auto'); % 'auto' or a number
compile       = ft_getopt(optarg, 'compile', 'no');   % can be 'auto', 'yes' or 'no'
backend       = ft_getopt(optarg, 'backend', []);     % the default will be determined by qsubfeval
queue         = ft_getopt(optarg, 'queue',   []);
submitoptions = ft_getopt(optarg, 'options', []);
batch         = ft_getopt(optarg, 'batch',   getbatch());               % this is a number that is automatically incremented
batchid       = ft_getopt(optarg, 'batchid', generatebatchid(batch));   % this is a string like user_host_pid_batch
display       = ft_getopt(optarg, 'display', 'no');
matlabcmd     = ft_getopt(optarg, 'matlabcmd', []);
jvm           = ft_getopt(optarg, 'jvm', 'yes');
whichfunction = ft_getopt(optarg, 'whichfunction');   % the complete filename to the function, including path
rerunable     = ft_getopt(optarg, 'rerunable');       % the default is determined in qsubfeval
sleep         = ft_getopt(optarg, 'sleep', 0.5);      % in seconds

% skip the optional key-value arguments
if ~isempty(optbeg)
  varargin = varargin(1:(optbeg-1));
end

if isstruct(fname)
  % the function has been compiled by qsubcompile
  fcomp = fname;
  % continue with the original function name
  fname = fcomp.fname;
else
  fcomp = [];
end

% determine which function it is
if isempty(whichfunction)
  if ischar(fname)
    whichfunction = which(fname);
  elseif isa(fname, 'function_handle')
    whichfunction = which(func2str(fname));
  end
end

% if the first attempt failed, it might be due a function that is private to the calling function
if isempty(whichfunction)
  s = dbstack('-completenames');
  s = s(2); % qsubcellfun is the first, the calling function is the second
  if ischar(fname)
    whichfunction = which(fullfile(fileparts(s.file), 'private', fname));
  elseif isa(fname, 'function_handle')
    whichfunction = which(fullfile(fileparts(s.file), 'private', func2str(fname)));
  end
  if ~isempty(whichfunction)
    warning('assuming %s as full function name', whichfunction);
  end
  clear s
end

% there are potentially errors to catch from the which() function
if isempty(whichfunction) && ischar(fname)
  error('not a valid M-file "%s"', fname);
end

% determine the number of input arguments and the number of jobs
numargin    = numel(varargin);
numjob      = numel(varargin{1});

% determine the number of MATLAB jobs to "stack" together into separate qsub jobs
if isequal(stack, 'auto')
  if ~isempty(timreq)
    stack = floor(180/timreq);
  else
    stack = 1;
  end
end

% ensure that the stacking is not higher than the number of jobs
stack = min(stack, numjob);

% give some feedback about the stacking
if stack>1
  fprintf('stacking %d MATLAB jobs in each qsub job\n', stack);
end

% prepare some arrays that are used for bookkeeping
jobid       = cell(1, numjob);
puttime     = nan(1, numjob);
timused     = nan(1, numjob);
memused     = nan(1, numjob);
submitted   = false(1, numjob);
collected   = false(1, numjob);
submittime  = inf(1, numjob);
collecttime = inf(1, numjob);

% it can be difficult to determine the number of output arguments
try
  if isequal(fname, 'cellfun') || isequal(fname, @cellfun)
    if isa(varargin{1}{1}, 'char') || isa(varargin{1}{1}, 'function_handle')
      numargout = nargout(varargin{1}{1});
    elseif isa(varargin{1}{1}, 'struct')
      % the function to be executed has been compiled
      fcomp = varargin{1}{1};
      numargout = nargout(fcomp.fname);
    end
  else
    numargout = nargout(fname);
  end
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

% running a compiled version in parallel takes no MATLAB licenses
% auto compilation will be attempted if the total batch takes more than 30 minutes
if (strcmp(compile, 'auto') && (numjob*timreq/3600)>0.5) || istrue(compile)
  try
    % try to compile into a stand-allone application
    fcomp = qsubcompile(fname, 'batch', batch, 'batchid', batchid);
  catch
    if istrue(compile)
      % the error that was caught is critical
      rethrow(lasterror);
    elseif strcmp(compile, 'auto')
      % compilation was only optional, the caught error is not critical
      warning(lasterr);
    end
  end % try-catch
end % if compile

if stack>1
  % combine multiple jobs in one, the idea is to use recursion like this
  % a = {{@plus, @plus}, {{1}, {2}}, {{3}, {4}}}
  % b = cellfun(@cellfun, a{:})

  % these options will be passed to the recursive call after being modified further down
  if ~any(strcmpi(optarg, 'timreq'))
    optarg{end+1} = 'timreq';
    optarg{end+1} = timreq;
  end
  if ~any(strcmpi(optarg, 'stack'))
    optarg{end+1} = 'stack';
    optarg{end+1} = stack;
  end
  if ~any(strcmpi(optarg, 'UniformOutput'))
    optarg{end+1} = 'UniformOutput';
    optarg{end+1} = UniformOutput;
  end
  if ~any(strcmpi(optarg, 'whichfunction'))
    optarg{end+1} = 'whichfunction';
    optarg{end+1} = whichfunction;
  end
  if ~any(strcmpi(optarg, 'compile'))
    optarg{end+1} = 'compile';
    optarg{end+1} = compile;
  end

  % update these settings for the recursive call
  optarg{find(strcmpi(optarg, 'timreq'))+1}        = timreq*stack;
  optarg{find(strcmpi(optarg, 'stack'))+1}         = 1;
  optarg{find(strcmpi(optarg, 'UniformOutput'))+1} = false;
  optarg{find(strcmpi(optarg, 'compile'))+1}       = false;

  % FIXME the partitioning can be further perfected
  partition     = floor((0:numjob-1)/stack)+1;
  numpartition  = partition(end);

  stackargin = cell(1,numargin+3); % include the fname, uniformoutput, false
  if istrue(compile)
    if ischar(fcomp.fname)
      % it should contain function handles, not strings
      stackargin{1} = repmat({str2func(fcomp.fname)}, 1, numpartition);
    else
      stackargin{1} = repmat({fcomp.fname}, 1, numpartition);
    end
  else
    if ischar(fname)
      % it should contain function handles, not strings
      stackargin{1} = repmat({str2func(fname)}, 1, numpartition);
    else
      stackargin{1} = repmat({fname}, 1, numpartition);
    end
  end
  stackargin{end-1} = repmat({'uniformoutput'},1,numpartition);  % uniformoutput
  stackargin{end}   = repmat({false},1,numpartition);            % false

  % reorganize the original input into the stacked format
  for i=1:numargin
    tmp = cell(1,numpartition);
    for j=1:numpartition
      tmp{j} = {varargin{i}{partition==j}};
    end
    stackargin{i+1} = tmp; % note that the first element is the fname
    clear tmp
  end

  stackargout = cell(1,numargout);
  [stackargout{:}] = qsubcellfun(@cellfun, stackargin{:}, optarg{:});

  % reorganise the stacked output into the original format
  for i=1:numargout
    tmp = cell(size(varargin{1}));
    for j=1:numpartition
      tmp(partition==j) = stackargout{i}{j};
    end
    varargout{i} = tmp;
    clear tmp
  end

  if numargout>0 && UniformOutput
    [varargout{:}] = makeuniform(varargout{:});
  end

  return;
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

for submit=1:numjob
  % redistribute the input arguments
  argin = cell(1, numargin);
  for j=1:numargin
    argin{j} = varargin{j}{submit};
  end

  % submit the job
  if ~isempty(fcomp)
    % use the compiled version
    [curjobid, curputtime] = qsubfeval(fcomp, argin{:}, 'memreq', memreq, 'timreq', timreq, 'memoverhead', memoverhead, 'timoverhead', timoverhead, 'diary', diary, 'batch', batch, 'batchid', batchid, 'backend', backend, 'options', submitoptions, 'queue', queue, 'matlabcmd', matlabcmd, 'display', display, 'jvm', jvm, 'nargout', numargout, 'whichfunction', whichfunction, 'rerunable', rerunable);
  else
    % use the non-compiled version
    [curjobid, curputtime] = qsubfeval(fname, argin{:}, 'memreq', memreq, 'timreq', timreq, 'memoverhead', memoverhead, 'timoverhead', timoverhead, 'diary', diary, 'batch', batch, 'batchid', batchid, 'backend', backend, 'options', submitoptions, 'queue', queue, 'matlabcmd', matlabcmd, 'display', display, 'jvm', jvm, 'nargout', numargout, 'whichfunction', whichfunction, 'rerunable', rerunable);
  end

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

      if isempty(argout) && StopOnError==false
        % this happens if an error was detected in qsubget and StopOnError is false
        % replace the output of the failed jobs with []
        argout = repmat({[]}, 1, numargout);
      end

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

  pausejava(sleep);
end % while

% ensure the output to have the same size/dimensions as the input
for i=1:numargout
  varargout{i} = reshape(varargout{i}, size(varargin{1}));
end

if numargout>0 && UniformOutput
  [varargout{:}] = makeuniform(varargout{:});
end

% clean up the remains of the compilation
if (strcmp(compile, 'yes') || strcmp(compile, 'auto')) && ~isempty(fcomp)
  % the extension might be .app or .exe or none
  system(sprintf('rm -rf %s',         fcomp.batchid)); % on Linux
  system(sprintf('rm -rf %s.app',     fcomp.batchid)); % on Apple OS X
  system(sprintf('rm -rf %s.exe',     fcomp.batchid)); % on Windows
  system(sprintf('rm -rf run_%s*.sh', fcomp.batchid));
end

% compare the time used inside this function with the total execution time
fprintf('computational time = %.1f sec, elapsed = %.1f sec, speedup %.1f x\n', nansum(timused), toc(stopwatch), nansum(timused)/toc(stopwatch));

if all(puttime>timused)
  warning('the job submission took more time than the actual execution');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = makeuniform(varargin)
varargout = varargin;
numargout = numel(varargin);
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
