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
%   backend     = string, can be 'torque', 'sge', 'slurm', 'lsf', 'system', 'local' (default is automatic)
%   diary       = string, can be 'always', 'never', 'warning', 'error' (default = 'error')
%   queue       = string, which queue to submit the job in (default is empty)
%   waitfor     = string or cell-array of strings, jobids of jobs to wait on finishing
%                 before executing the current job (default is empty)
%   options     = string, additional options that will be passed to qsub/srun (default is empty)
%   batch       = number, of the bach to which the job belongs. When called by QSUBCELLFUN
%                 it will be a number that is automatically incremented over subsequent calls.
%   batchid     = string that is used for the compiled application filename and to identify
%                 the jobs in the queue, the default is automatically determined and looks
%                 like user_host_pid_batch.
%   matlabcmd   = string, the Linux command line to start MATLAB on the compute nodes (default is automatic
%   display     = 'yes' or 'no', whether the nodisplay option should be passed to MATLAB (default = 'no', meaning nodisplay)
%   jvm         = 'yes' or 'no', whether the nojvm option should be passed to MATLAB (default = 'yes', meaning with jvm)
%   rerunable   = 'yes' or 'no', whether the job can be restarted on a torque/maui/moab cluster (default = 'no')
%
% See also QSUBCELLFUN, QSUBGET, FEVAL, DFEVAL, DFEVALASYNC

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

% these are used to speed up the processing of multiple function calls with
% the same input arguments (e.g. from peercellfun)
persistent previous_argin
persistent previous_matlabcmd

% keep track of the time
stopwatch = tic;

% convert the input arguments into something that strcmp can work with
strargin = varargin;
strargin(~cellfun(@ischar, strargin)) = {''};

% locate the begin of the optional key-value arguments
optbeg = false(size(strargin));
optbeg = optbeg | strcmp('memreq',        strargin);
optbeg = optbeg | strcmp('timreq',        strargin);
optbeg = optbeg | strcmp('diary',         strargin);
optbeg = optbeg | strcmp('batch',         strargin);
optbeg = optbeg | strcmp('batchid',       strargin);
optbeg = optbeg | strcmp('timoverhead',   strargin);
optbeg = optbeg | strcmp('memoverhead',   strargin);
optbeg = optbeg | strcmp('backend',       strargin);
optbeg = optbeg | strcmp('queue',         strargin);
optbeg = optbeg | strcmp('options',       strargin);
optbeg = optbeg | strcmp('matlabcmd',     strargin);
optbeg = optbeg | strcmp('jvm',           strargin);
optbeg = optbeg | strcmp('display',       strargin);
optbeg = optbeg | strcmp('nargout',       strargin);
optbeg = optbeg | strcmp('whichfunction', strargin);
optbeg = optbeg | strcmp('waitfor',       strargin);
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
backend       = ft_getopt(optarg, 'backend', []);                 % the defaultbackend helper function will be used to determine the default
queue         = ft_getopt(optarg, 'queue', []);                   % the default is specified further down in the code
submitoptions = ft_getopt(optarg, 'options', []);
display       = ft_getopt(optarg, 'display', 'no');
matlabcmd     = ft_getopt(optarg, 'matlabcmd', []);
jvm           = ft_getopt(optarg, 'jvm', 'yes');
numargout     = ft_getopt(optarg, 'nargout', []);
whichfunction = ft_getopt(optarg, 'whichfunction');               % the complete filename to the function, including path
waitfor       = ft_getopt(optarg, 'waitfor', {});                 % default is empty cell-array
rerunable     = ft_getopt(optarg, 'rerunable');                   % the default is determined in qsubfeval

% skip the optional key-value arguments
if ~isempty(optbeg)
  varargin = varargin(1:(optbeg-1));
end

if isempty(backend)
  % use the system default backend
  backend = defaultbackend;
else
  % use the user-specified backend
  % this makes it persistent and available to qsublist
  defaultbackend(backend);
end

% it should be specified as a cell-array of strings
if ischar(waitfor)
  waitfor = {waitfor};
end
% remove empty elements
waitfor = waitfor(~cellfun(@isempty, waitfor));

% determine whether the function has been compiled
compiled = isstruct(varargin{1});

if compiled
  % the function has been compited by qsubcompile
  compiledfun = varargin{1}.executable;
  % continue with the original function name
  varargin{1} = varargin{1}.fname;
end

hostname = gethostname();
if isempty(queue) && ~compiled && (~isempty(regexp(hostname, '^dccn-c', 'once')) || ~isempty(regexp(hostname, '^mentat', 'once')))
  % At the DCCN we want the non-compiled distributed MATLAB jobs to be queued in the "matlab" queue. This
  % routes them to specific multi-core machines and limits the number of licenses that can be claimed at once.
  queue = 'matlab';
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
options = {'pwd', curPwd, 'path', getcustompath, 'global', getglobal, 'diary', diary, 'memreq', memreq, 'timreq', timreq, 'randomseed', randomseed, 'nargout', numargout, 'whichfunction', whichfunction, 'rerunable', rerunable};

inputfile    = fullfile(curPwd, sprintf('%s_input.mat', jobid));
matlabscript = fullfile(curPwd, sprintf('%s.m', jobid));

% rename and save the variables
argin = varargin;
optin = options;
s1 = whos('argin');
s2 = whos('optin');
% if variables < ~1 GB, store it in old (uncompressed) format, which is faster
if (s1.bytes + s2.bytes < 1024^3)
  save(inputfile, 'argin', 'optin', '-v6');
else
  save(inputfile, 'argin', 'optin', '-v7.3');
end

if ~compiled
  
  if ~isempty(matlabcmd)
    % take the user-specified matlab startup script
  elseif isempty(previous_matlabcmd)
    % determine the name of the matlab startup script
    
    if ft_platform_supports('program_invocation_name')
      % supported in GNU Octave
      matlabcmd = program_invocation_name();
    elseif ~isempty(getenv('MATLAB_BIN'))
      % supported on Linux + R2014b, perhaps also on others
      matlabcmd = getenv('MATLAB_BIN');
    elseif ~isempty(getenv('MATLABDIR'))
      % supported on Linux + R2012b, perhaps also on others
      matlabcmd = fullfile(getenv('MATLABDIR'), 'bin/matlab');
    else
      matlabcmd = '';
      % try all versions between 7.1 and 7.9
      for matlab_version=71:79
        matlab_version_decimated=matlab_version*.1;
        if ft_platform_supports('matlabversion',matlab_version_decimated)
          matlabcmd = sprintf('matlab%d',matlab_version);
          break;
        end
      end
      if isempty(matlabcmd)
        matlabcmd = sprintf('matlab%s', version('-release')); % the version command returns a string like '2014a'
      end
    end
    
    if system(sprintf('which %s > /dev/null', matlabcmd))==1
      % the linux command "which" returns 0 on succes and 1 on failure
      warning('the executable for "%s" could not be found, trying "matlab" instead', matlabcmd);
      % use whatever is available as default
      matlabcmd = 'matlab';
    end
    
    % keep the matlab command for subsequent calls, this will
    % avoid subsequent attempts to set the matlabcmd
    % and the system('which ...') call on the scheduling of subsequent
    % distributed jobs
    previous_matlabcmd = matlabcmd;
  else
    % re-use the matlab command that was determined on the previous call to this function
    matlabcmd = previous_matlabcmd;
  end
  
  if ft_platform_supports('singleCompThread')
    % this is only supported for version 7.8 onward
    matlabcmd = [matlabcmd ' -singleCompThread'];
  end
  
  % these options can be appended regardless of the version
  if ft_platform_supports('nosplash');
    matlabcmd = [matlabcmd ' -nosplash'];
  end
  if ~istrue(display)
    if ft_platform_supports('nodisplay');
      % Matlab
      matlabcmd = [matlabcmd ' -nodisplay'];
    end
    if ft_platform_supports('no-gui');
      % GNU Octave
      matlabcmd = [matlabcmd ' --no-gui'];
    end
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
  
end % if ~compiled

% set the job requirements according to the users specification
switch backend
  case 'local'
    % this is for testing the execution in case no cluster is available,
    % for example when working on the road with a laptop
    
    cmdline = [];
    
  case 'system'
    % this is for testing the execution in case no cluster is available,
    % for example when working on the road with a laptop
    
    if compiled
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
      submitoptions = [submitoptions sprintf('-l h_rt=%.0f ', timreq+timoverhead)];
    end
    
    if ~isempty(memreq) && ~isnan(memreq) && ~isinf(memreq)
      submitoptions = [submitoptions sprintf('-l mem_free=%.0fG ', round((memreq+memoverhead)/1024^3))];
    end
    
    if compiled
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
      submitoptions = [submitoptions sprintf(' -l walltime=%.0f ', timreq+timoverhead)];
    end
    
    if ~isempty(memreq) && ~isnan(memreq) && ~isinf(memreq)
      % mem is the real memory, vmem is the virtual, pmem and pvmem relate to the memory per process in case of an MPI job with multiple processes
      submitoptions = [submitoptions sprintf(' -l mem=%.0f ',   memreq+memoverhead)];
      %   submitoptions = [submitoptions sprintf(' -l vmem=%.0f ',  memreq+memoverhead)];
      %   submitoptions = [submitoptions sprintf(' -l pmem=%.0f ',  memreq+memoverhead)];
      %   submitoptions = [submitoptions sprintf(' -l pvmem=%.0f ', memreq+memoverhead)];
    end
    
    if ~isempty(waitfor)
      % waitfor contains the jobids of the jobs to wait for
      submitoptions = [submitoptions '-W depend=afterok'];
      for iJob = 1:numel(waitfor)
        submitoptions = [submitoptions sprintf(':%s',qsublist('getpbsid', waitfor{iJob}))];
      end
    end
    
    % In the command below both stderr and stout are redirected to /dev/null,
    % so any output information will not be available for inspection.
    % However, any matlab errors will be reported back by fexec.
    % cmdline = ['qsub -e /dev/null -o /dev/null -N ' jobid ' ' requirements shellscript];
    
    if compiled
      % create the command line for the compiled application
      cmdline = sprintf('%s %s %s', compiledfun, matlabroot, jobid);
    else
      % create the shell commands to execute matlab
      cmdline = sprintf('%s -r \\"%s\\"', matlabcmd, matlabscript);
    end
    
    if any(curPwd==' ')
      % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1898
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
    
    if compiled
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
    
    
  case 'lsf'
    % this is for Platform Load Sharing Facility (LSF)
    
    if isempty(submitoptions)
      % start with an empty string
      submitoptions = '';
    end
    
    if ~isempty(queue)
      submitoptions = [submitoptions sprintf('-q %s ', queue)];
    end
    
    if ~isempty(timreq) && ~isnan(timreq) && ~isinf(timreq)
      submitoptions = [submitoptions sprintf('-W %.0f ', ceil((timreq+timoverhead) / 60))]; % in minutes
    end
    
    if ~isempty(memreq) && ~isnan(memreq) && ~isinf(memreq)
      submitoptions = [submitoptions sprintf('-M %.0f ', ceil((memreq+memoverhead) / 1024^2))];  % in MB
    end
    
    % specifying the o and e names might be useful for the others as well
    logout = fullfile(curPwd, sprintf('%s.o', jobid));
    logerr = fullfile(curPwd, sprintf('%s.e', jobid));
    
    if compiled
      % create the command line for the compiled application
      cmdline = sprintf('%s %s %s', compiledfun, matlabroot, jobid);
    else
      % create the shell commands to execute matlab
      cmdline = sprintf('%s -r \\"%s\\"', matlabcmd, matlabscript);
    end
    
    % pass the command to qsub with all requirements
    cmdline = sprintf('echo "%s" | bsub -J %s %s -o %s -e %s', cmdline, jobid, submitoptions, logout, logerr);
    
  otherwise
    error('unsupported backend "%s"', backend);
    
end % switch

fprintf('submitting job %s...', jobid); % note the lack of the end-of-line, the qsub output will follow

if ~strcmp(backend, 'local')
  % the system call will also print some information to screen to complete the line
  [status, result] = system(cmdline);
  if status
    % this should have returned 0, the screen output in result will probably be informative
    error(result);
  end
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
  case 'lsf'
    % the result of bsub returns a string in format: "Job <job_number> is submitted to default queue <queue_name>"
    % parse the job number
    pbsid_beg = strfind(result, '<');
    pbsid_end = strfind(result, '>');
    result = result(pbsid_beg(1)+1:pbsid_end(1)-1);
  case 'sge'
    % in sge, the return string is "Your job <job_number> (<job_name>) has been submitted"
    result_words = tokenize(result, ' ');
    result = result_words{3};
  otherwise
    % for torque, it is enough to remove the white space
    result = strtrim(result);
end

fprintf(' %s id %s\n', backend, result);

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
