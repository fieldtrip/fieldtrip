function [jobid, puttime] = qsubfeval(varargin)

% QSUBFEVAL evaluates the specified MATLAB function on the input arguments
% using the Torque, SGE, PBS or SLURM batch queue system.
%
% Use as
%   jobid  = qsubfeval(fname, arg1, arg2, ...)
%   argout = qsubget(jobid, ...)
%
% This function has a number of optional arguments that have to passed
% as key-value pairs at the end of the list of input arguments. All other
% input arguments (including other key-value pairs) will be passed to the
% function to be evaluated.
%   memreq      = number in bytes, how much memory does the job require (no default)
%   memoverhead = number in bytes, how much memory to account for MATLAB itself (default = 1024^3, i.e. 1GB)
%   timreq      = number in seconds, how much time does the job require (no default)
%   timoverhead = number in seconds, how much time to allow MATLAB to start (default = 180 seconds)
%   backend     = string, can be 'sge', 'torque', 'slurm', 'system', 'local' (default is automatic)
%   diary       = string, can be 'always', 'never', 'warning', 'error' (default = 'error')
%   queue       = string, which queue to submit the job in (default is empty)
%   options     = string, additional options that will be passed to qsub/srun (default is empty)
%   batch       = number, of the bach to which the job belongs. When called by QSUBCELLFUN
%                 it will be a number that is automatically incremented over subsequent calls.
%   batchid     = string that is used for the compiled application filename and to identify
%                 the jobs in the queue, the default is automatically determined and looks
%                 like user_host_pid_batch.
%   display     = 'yes' or 'no', whether the nodisplay option should be passed to MATLAB (default = 'no', meaning nodisplay)
%   jvm         = 'yes' or 'no', whether the nojvm option should be passed to MATLAB (default = 'yes', meaning with jvm)
%
% See also QSUBCELLFUN, QSUBGET, FEVAL, DFEVAL, DFEVALASYNC

% -----------------------------------------------------------------------
% Copyright (C) 2011-2012, Robert Oostenveld
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

% these are used to speed up the processing of multiple function calls with
% the same input arguments (e.g. from peercellfun)
persistent previous_argin
persistent previous_matlabcmd

% keep track of the time
stopwatch = tic;

% check if torque or sge is present and running
if ~isempty(getenv('SGE_ROOT'))
  defaultbackend = 'sge';
elseif ~isempty(getenv('TORQUEHOME'))
  defaultbackend = 'torque';
elseif ~isempty(getenv('CONDOR_ARCH'))
  % this has not been tested and I am not 100% sure that this is the right variable to probe
  defaultbackend = 'condor';
elseif ~isempty(getenv('SLURM_ENABLE'))
  defaultbackend = 'slurm';
else
  % backend=local causes the job to be executed in this MATLAB by feval
  % backend=system causes the job to be executed on the same computer using system('matlab -r ...')
  defaultbackend = 'local';
end

hostname = gethostname();
if ~isempty(regexp(hostname, '^dccn-c', 'once')) || ~isempty(regexp(hostname, '^mentat', 'once'))
  % At the DCCN we want the distributed MATLAB jobs to be queued in the
  % "matlab" queue. This routes them to specific multi-core machines and
  % limits the number of licenses that can be claimed at once.
  defaultqueue = 'matlab';
else
  % let the queueing system decide
  defaultqueue = [];
end

% convert the input arguments into something that strmatch can work with
strargin = varargin;
strargin(~cellfun(@ischar, strargin)) = {''};

% locate the begin of the optional key-value arguments
optbeg = false(size(strargin));
optbeg = optbeg | strcmp('memreq',      strargin);
optbeg = optbeg | strcmp('timreq',      strargin);
optbeg = optbeg | strcmp('diary',       strargin);
optbeg = optbeg | strcmp('batch',       strargin);
optbeg = optbeg | strcmp('batchid',     strargin);
optbeg = optbeg | strcmp('timoverhead', strargin);
optbeg = optbeg | strcmp('memoverhead', strargin);
optbeg = optbeg | strcmp('backend',     strargin);
optbeg = optbeg | strcmp('queue',       strargin);
optbeg = optbeg | strcmp('options',     strargin);
optbeg = optbeg | strcmp('jvm',         strargin);
optbeg = optbeg | strcmp('display',     strargin);
optbeg = optbeg | strcmp('nargout',     strargin);
optbeg = find(optbeg);
optarg = varargin(optbeg:end);

% check the required input arguments
ft_checkopt(optarg, 'memreq', 'numericscalar');
ft_checkopt(optarg, 'timreq', 'numericscalar');

% get the optional input arguments
memreq        = ft_getopt(optarg, 'memreq');
timreq        = ft_getopt(optarg, 'timreq');
diary         = ft_getopt(optarg, 'diary');
batch         = ft_getopt(optarg, 'batch', 1);
batchid       = ft_getopt(optarg, 'batchid');
timoverhead   = ft_getopt(optarg, 'timoverhead', 180);            % allow some overhead to start up the MATLAB executable
memoverhead   = ft_getopt(optarg, 'memoverhead', 1024*1024*1024); % allow some overhead for the MATLAB executable in memory
backend       = ft_getopt(optarg, 'backend', defaultbackend);     % can be torque, local, sge
queue         = ft_getopt(optarg, 'queue', defaultqueue);
submitoptions = ft_getopt(optarg, 'options', []);
display       = ft_getopt(optarg, 'display', 'no');
jvm           = ft_getopt(optarg, 'jvm', 'yes');
numargout     = ft_getopt(optarg, 'nargout', []);

% skip the optional key-value arguments
if ~isempty(optbeg)
  varargin = varargin(1:(optbeg-1));
end

% determine whether the function has been compiled
compile = isstruct(varargin{1});

if isa(varargin{1}, 'struct')
  % the function has been compited by qsubcompile
  compiledfun = varargin{1}.executable;
  % continue with the original function name
  varargin{1} = varargin{1}.fname;
end

if ~isempty(previous_argin) && ~isequal(varargin{1}, previous_argin{1})
  % this can be skipped if the previous call used the same function
  if ischar(varargin{1}) && isempty(which(varargin{1}))
    error('Not a valid M-file (%s).', varargin{1});
  end
end

% create a unique identifier for the job (string)
jobid = generatejobid(batch, batchid);

% get the current working directory to store the temp files in
curPwd = getcustompwd();

% each job should have a different random number sequence
randomseed = rand(1)*double(intmax);

% pass some options that influence the remote execution
options = {'pwd', curPwd, 'path', getcustompath, 'global', getglobal, 'diary', diary, 'memreq', memreq, 'timreq', timreq, 'randomseed', randomseed, 'nargout', numargout};

inputfile    = fullfile(curPwd, sprintf('%s_input.mat', jobid));
matlabscript = fullfile(curPwd, sprintf('%s.m', jobid));

% rename and save the variables
argin = varargin;
optin = options;
save(inputfile, 'argin', 'optin');

if ~compile
  
  if isempty(previous_matlabcmd)
    % determine the name of the matlab startup script
    if matlabversion(7.1)
      matlabcmd = 'matlab71';
    elseif matlabversion(7.2)
      matlabcmd = 'matlab72';
    elseif matlabversion(7.3)
      matlabcmd = 'matlab73';
    elseif matlabversion(7.4)
      matlabcmd = 'matlab74';
    elseif matlabversion(7.5)
      matlabcmd = 'matlab75';
    elseif matlabversion(7.6)
      matlabcmd = 'matlab76';
    elseif matlabversion(7.7)
      matlabcmd = 'matlab77';
    elseif matlabversion(7.8) % 2009a
      matlabcmd = 'matlab78';
    elseif matlabversion(7.9) % 2009b
      matlabcmd = 'matlab79';
    elseif matlabversion('2010a')
      matlabcmd = 'matlab2010a';
    elseif matlabversion('2010b')
      matlabcmd = 'matlab2010b';
    elseif matlabversion('2011a')
      matlabcmd = 'matlab2011a';
    elseif matlabversion('2011b')
      matlabcmd = 'matlab2011b';
    elseif matlabversion('2012a')
      matlabcmd = 'matlab2012a';
    elseif matlabversion('2012b')
      matlabcmd = 'matlab2012b';
    else
      % use whatever is available as default
      matlabcmd = 'matlab';
    end
    
    if system(sprintf('which %s > /dev/null', matlabcmd))==1
      % the linux command "which" returns 0 on succes and 1 on failure
      warning('the executable for "%s" could not be found, trying "matlab" instead', matlabcmd);
      % use whatever is available as default
      matlabcmd = 'matlab';
    end
    
    % keep the matlab command for subsequent calls, this will save all the matlabversion calls
    % and the system('which ...') call on the scheduling of subsequent distributed jobs
    previous_matlabcmd = matlabcmd;
  else
    % re-use the matlab command that was determined on the previous call to this function
    matlabcmd = previous_matlabcmd;
  end
  
  if matlabversion(7.8, inf)
    % this is only supported for version 7.8 onward
    matlabcmd = [matlabcmd ' -singleCompThread'];
  end
  
  % these options can be appended regardless of the version
  matlabcmd = [matlabcmd ' -nosplash'];
  if ~istrue(display)
    matlabcmd = [matlabcmd ' -nodisplay'];
  end
  if ~istrue(jvm)
    matlabcmd = [matlabcmd ' -nojvm'];
  end
  
  % create the matlab script commands (one entry per line)
  matlabscript = [...
    'restoredefaultpath;',...
    sprintf('addpath(''%s'');', fileparts(mfilename('fullpath'))),...
    sprintf('qsubexec(''%s'');', fullfile(pwd, jobid)),...
    sprintf('exit')];

end % if ~compile

% set the job requirements according to the users specification
switch backend
  case 'local'
    % this is for testing the execution in case no cluster is available,
    % for example when working on the road with a laptop
    
    cmdline = [];
    
  case 'system'
    % this is for testing the execution in case no cluster is available,
    % for example when working on the road with a laptop
    
    if compile
      % create the command line for the compiled application
      cmdline = sprintf('%s %s %s', compiledfun, matlabroot, jobid);
    else
      % create the shell commands to execute matlab
      cmdline = sprintf('%s -r "%s"', matlabcmd, matlabscript);
    end
    
  case 'sge'
    % this is for Sun Grid Engine, Oracle Grid Engine, and other derivatives
    
    if isempty(submitoptions)
      % start with an empty string
      submitoptions = '';
    end
    
    if ~isempty(queue)
      submitoptions = [submitoptions sprintf('-q %s ', queue)];
    end
    
    if ~isempty(timreq) && ~isnan(timreq) && ~isinf(timreq)
      submitoptions = [submitoptions sprintf('-l h_rt=%d ', timreq+timoverhead)];
    end
    
    if ~isempty(memreq) && ~isnan(memreq) && ~isinf(memreq)
      submitoptions = [submitoptions sprintf('-l h_vmem=%.0f ', memreq+memoverhead)];
    end
    
    if compile
      % create the command line for the compiled application
      cmdline = sprintf('%s %s %s', compiledfun, matlabroot, jobid);
    else
      % create the shell commands to execute matlab
      cmdline = sprintf('%s -r \\"%s\\"', matlabcmd, matlabscript);
    end
    
    % pass the command to qsub with all requirements
    cmdline = sprintf('echo "%s" | qsub -N %s %s -cwd -o %s -e %s', cmdline, jobid, submitoptions, curPwd, curPwd);
    
  case 'torque'
    % this is for PBS, Torque, and other derivatives
    
    if isempty(submitoptions)
      % start with an empty string
      submitoptions = '';
    end
    
    if ~isempty(queue)
      submitoptions = [submitoptions sprintf(' -q %s ', queue)];
    end
    
    if ~isempty(timreq) && ~isnan(timreq) && ~isinf(timreq)
      submitoptions = [submitoptions sprintf(' -l walltime=%d ', timreq+timoverhead)];
    end
    
    if ~isempty(memreq) && ~isnan(memreq) && ~isinf(memreq)
      % mem is the real memory, vmem is the virtual, pmem and pvmem relate to the memory per process in case of an MPI job with multiple processes
      submitoptions = [submitoptions sprintf(' -l mem=%.0f ',   memreq+memoverhead)];
      %   submitoptions = [submitoptions sprintf(' -l vmem=%.0f ',  memreq+memoverhead)];
      %   submitoptions = [submitoptions sprintf(' -l pmem=%.0f ',  memreq+memoverhead)];
      %   submitoptions = [submitoptions sprintf(' -l pvmem=%.0f ', memreq+memoverhead)];
    end
    
    % In the command below both stderr and stout are redirected to /dev/null,
    % so any output information will not be available for inspection.
    % However, any matlab errors will be reported back by fexec.
    % cmdline = ['qsub -e /dev/null -o /dev/null -N ' jobid ' ' requirements shellscript];
    
    if compile
      % create the command line for the compiled application
      cmdline = sprintf('%s %s %s', compiledfun, matlabroot, jobid);
    else
      % create the shell commands to execute matlab
      cmdline = sprintf('%s -r \\"%s\\"', matlabcmd, matlabscript);
    end
    
    if any(curPwd==' ')
      % see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1898
      error('you cannot execute jobs from within a directory that has a space in its name');
    end
    
    % pass the command to qsub with all requirements
    cmdline = sprintf('echo "%s" | qsub -N %s %s -d "%s" -o "%s" -e "%s"', cmdline, jobid, submitoptions, curPwd, curPwd, curPwd);
    
  case 'slurm'
    % this is for Simple Linux Utility for Resource Management
    
    if isempty(submitoptions)
      % start with an empty string
      submitoptions = '';
    end
    
    if ~isempty(queue)
      % with slurm queues are "partitions"
      submitoptions = [submitoptions sprintf(' --partition=%s ', queue)];
    end
    
    if ~isempty(timreq) && ~isnan(timreq) && ~isinf(timreq)
      % TESTME this is experimental and needs more testing!
      % submitoptions = [submitoptions sprintf('--time=%d ', timreq+timoverhead)];
    end
    
    if ~isempty(memreq) && ~isnan(memreq) && ~isinf(memreq)
      % TESTME this is experimental and needs more testing!
      % submitoptions = [submitoptions sprintf('--mem-per-cpu=%.0f ', round((memreq+memoverhead)./1024^2))];
    end
    
    % specifying the o and e names might be useful for the others as well
    logout = fullfile(curPwd, sprintf('%s.o', jobid));
    logerr = fullfile(curPwd, sprintf('%s.e', jobid));
    
    if compile
      % create the command line for the compiled application
      cmdline = sprintf('%s %s %s', compiledfun, matlabroot, jobid);
    else
      % create the shell commands to execute matlab
      % we decided to use srun instead of sbatch since handling job paramters is easier this way
      %
      % nohup was found to signficantly speed up the submission. Due to the existing error handling its safe to detach to the init, but debugging
      % gets harder since output will be redirected to nohpu.out and thus overwritten everytime qsubfeval is launched. Using nohup only makes sense
      % if you intend to sumbit jobs which compute in less than a minute since the difference in submit time is about 3-4 seconds per job only!
      % cmdline = sprintf('nohup srun --job-name=%s %s --output=%s --error=%s %s -r "%s" & ', jobid, submitoptions, logout, logerr, matlabcmd, matlabscript);
      cmdline = sprintf('srun --job-name=%s %s --output=%s --error=%s %s -r "%s" ', jobid, submitoptions, logout, logerr, matlabcmd, matlabscript);
    end
    
  case 'condor'
    % this is highly experimental and contains some first ideas following the discussion with Rhodri
    
    % create a condor submit script
    submitfile = fullfile(curPwd, sprintf('%s.condor', jobid));
    
    % the Condor submit script should look something like this
    fid = fopen(submitfile, 'wt');
    fprintf(fid, '# Condor submit script\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Executable     = %s\n', matlabcmd);
    fprintf(fid, 'Arguments      = -r "%s"\n', matlabscript);
    % the timreq and memrequ should be inserted here
    fprintf(fid, 'Requirements   = Memory >= 32 && OpSys == "LINUX" && Arch =="INTEL"\n');
    fprintf(fid, 'Rank           = Memory >= 64\n');
    fprintf(fid, 'Image_Size     = 28 Meg\n');
    fprintf(fid, '\n');
    % these output files should match with the ones expected in qsubget
    fprintf(fid, 'Error   = %s.err\n', jobid);
    fprintf(fid, 'Output  = %s.out\n', jobid);
    fprintf(fid, 'Log     = %s.log\n', jobid);
    fprintf(fid, '\n');
    fprintf(fid, 'Queue\n');
    fclose(fid);
    
    cmdline = sprintf('condor_submit %s', submitfile);
    
  otherwise
    error('unsupported backend "%s"', backend);
    
end % switch

fprintf('submitting job %s...', jobid); % note the lack of the end-of-line, the qsub outpt will follow
if ~strcmp(backend, 'local')
  [status, result] = system(cmdline);
else
  % this will read the job input *.mat file, call feval with all try-catch
  % precautions, measure time and memory and eventually write the results to
  % the job output *.mat file
  fprintf('\n');
  qsubexec(jobid);
end

switch backend
  case 'slurm'
    % srun will not return a jobid (besides in verbose mode) we decided to use jobname=jobid instead to identify processes
    % since the jobid is a uniq identifier for every job!
    result = jobid;
  case 'local'
    % the job was executed by a local feval call, but the results will still be written in a job file
    result = jobid;
  otherwise
    % for torque and sge it is enough to remove the white space
    result = strtrim(result);
end

fprintf(' qstat job id %s\n', result);

% both Torque and SGE will return a log file with stdout and stderr information
% for local execution we have to emulate these files, because qsubget expects them
if any(strcmp(backend, {'system', 'local'}))
  logout       = fullfile(curPwd, sprintf('%s.o', jobid));
  logerr       = fullfile(curPwd, sprintf('%s.e', jobid));
  fclose(fopen(logout, 'w'));
  fclose(fopen(logerr, 'w'));
end

% add the job to the persistent list, this is used for cleanup in case of Ctrl-C
qsublist('add', jobid, result);

puttime = toc(stopwatch);

% remember the input arguments to speed up subsequent calls
previous_argin  = varargin;

