function [jobid, puttime] = qsubfeval(varargin)

% QSUBFEVAL evaluates the specified MATLAB function on the input arguments
% using the Torque or SGE batch queue system.
%
% Use as
%   jobid  = qsubfeval(fname, arg1, arg2, ...)
%   argout = qsubget(jobid, ...)
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
% -----------------------------------------------------------------------

% these are used to speed up the processing of multiple function calls with
% the same input arguments (e.g. from peercellfun)
persistent previous_argin

% keep track of the time
stopwatch = tic;

% check if torque or sge is present and running
if ~isempty(getenv('SGE_ROOT'))
  defaultbackend = 'sge';
elseif ~isempty(getenv('TORQUEHOME'))
  defaultbackend = 'torque';
else
  defaultbackend = 'local';
end

% convert the input arguments into something that strmatch can work with
strargin = varargin;
strargin(~cellfun(@ischar, strargin)) = {''};

% locate the begin of the optional key-value arguments
optbeg = false(size(strargin));
optbeg = optbeg | strcmp('memreq',  strargin);
optbeg = optbeg | strcmp('timreq',  strargin);
optbeg = optbeg | strcmp('diary',   strargin);
optbeg = optbeg | strcmp('batch',   strargin);
optbeg = optbeg | strcmp('timoverhead', strargin);
optbeg = optbeg | strcmp('memoverhead', strargin);
optbeg = optbeg | strcmp('backend', strargin);
optbeg = find(optbeg);
optarg = varargin(optbeg:end);

% check the required input arguments
ft_checkopt(optarg, 'memreq', 'numericscalar');
ft_checkopt(optarg, 'timreq', 'numericscalar');

% get the optional input arguments
memreq      = ft_getopt(optarg, 'memreq');
timreq      = ft_getopt(optarg, 'timreq');
diary       = ft_getopt(optarg, 'diary',   []);
batch       = ft_getopt(optarg, 'batch',    1);
timoverhead = ft_getopt(optarg, 'timoverhead', 180);            % allow some overhead to start up the MATLAB executable
memoverhead = ft_getopt(optarg, 'memoverhead', 1024*1024*1024); % allow some overhead for the MATLAB executable in memory
backend     = ft_getopt(optarg, 'backend', defaultbackend);     % can be torque, local, sge

% skip the optional key-value arguments
if ~isempty(optbeg)
  varargin = varargin(1:(optbeg-1));
end

% determine whether the function has been compiled
compile = isstruct(varargin{1});

if isa(varargin{1}, 'function_handle')
  % convert the function handle back into a string (e.g. @plus should be 'plus')
  varargin{1} = func2str(varargin{1});
elseif isa(varargin{1}, 'struct')
  % the function has been compited by qsubcompile
  compiledfun = varargin{1}.executable;
  % continue with the original function name
  varargin{1} = varargin{1}.fname;
end

if ~isempty(previous_argin) && ~isequal(varargin{1}, previous_argin{1})
  % this can be skipped if the previous call used the same function
  if isempty(which(varargin{1}))
    error('Not a valid M-file (%s).', varargin{1});
  end
end

% create a unique identifier for the job (string)
jobid = generatejobid(batch);

% get the current working directory to store the temp files in
curPwd = getcustompwd();

% each job should have a different random number sequence
randomseed = rand(1)*double(intmax);

% pass some options that influence the remote execution
options = {'pwd', curPwd, 'path', getcustompath, 'global', getglobal, 'diary', diary, 'memreq', memreq, 'timreq', timreq, 'randomseed', randomseed};

inputfile    = fullfile(curPwd, sprintf('%s_input.mat', jobid));
shellscript  = fullfile(curPwd, sprintf('%s.sh', jobid));
matlabscript = fullfile(curPwd, sprintf('%s.m', jobid));

% rename and save the variables
argin = varargin;
optin = options;
save(inputfile, 'argin', 'optin');

if ~compile
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
  
  if system(sprintf('which %s', matlabcmd))==1
    % the linux command "which" returns 0 on succes and 1 on failure
    warning('the executable for "%s" could not be found, trying "matlab" instead', matlabcmd);
    % use whatever is available as default
    matlabcmd = 'matlab';
  end
  
  if matlabversion(7.8, inf)
    % this is only supported for version 7.8 onward
    matlabcmd = [matlabcmd ' -singleCompThread'];
  end
  
  % these options can be appended regardless of the version
  matlabcmd = [matlabcmd ' -nosplash -nodisplay'];
  
  % create the matlab script commands (one entry per line)
  matlabscript = [...
    'restoredefaultpath;',...
    sprintf('cd(''%s'');', curPwd),...
    sprintf('addpath(''%s'');', fileparts(mfilename('fullpath'))),...
    sprintf('qsubexec(''%s'');', jobid),...
    sprintf('exit')];
  
end % if ~compile

% set the job requirements according to the users specification
switch backend
  case 'local'
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
    requirements = '';
    if ~isempty(timreq)
      requirements = [requirements sprintf('-l h_rt=%d ', timreq+timoverhead)];
    end
    if ~isempty(memreq)
      requirements = [requirements sprintf('-l h_vmem=%.0f ',   memreq+memoverhead)];
    end
    
    if compile
      % create the command line for the compiled application
      cmdline = sprintf('%s %s %s', compiledfun, matlabroot, jobid);
    else
      % create the shell commands to execute matlab
      cmdline = sprintf('%s -r \\"%s\\"', matlabcmd, matlabscript);
    end
    
    % pass the command to qsub with all requirements
    cmdline = sprintf('echo "%s" | qsub -N %s %s -cwd -o %s -e %s', cmdline, jobid, requirements, curPwd, curPwd);
    
  case 'torque'
    % this is for PBS, Torque, and other derivatives
    requirements = '';
    if ~isempty(timreq)
      requirements = [requirements sprintf('-l walltime=%d ', timreq+timoverhead)];
    end
    if ~isempty(memreq)
      % mem is the real memory, vmem is the virtual, pmem and pvmem relate to the memory per process in case of an MPI job with multiple processes
      requirements = [requirements sprintf('-l mem=%.0f ',   memreq+memoverhead)];
      %   requirements = [requirements sprintf('-l vmem=%.0f ',  memreq+memoverhead)];
      %   requirements = [requirements sprintf('-l pmem=%.0f ',  memreq+memoverhead)];
      %   requirements = [requirements sprintf('-l pvmem=%.0f ', memreq+memoverhead)];
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
    
    % pass the command to qsub with all requirements
    cmdline = sprintf('echo "%s" | qsub -N %s %s -d %s -o %s -e %s', cmdline, jobid, requirements, curPwd, curPwd, curPwd);
    
  case 'slurm'
    % this is for Simple Linux Utility for Resource Management
    error('not yet implemented');
    
  otherwise
    error('unsupported backend "%s"', backend);
    
end % switch

fprintf('submitting job %s...', jobid);
[status, result] = system(cmdline);
fprintf(' qstat job id %s\n', strtrim(result));

% both Torque and SGE will return a log file with stdout and stderr information
% for local execution we have to emulate these files, because qsubget expects them
if strcmp(backend, 'local')
  logout       = fullfile(curPwd, sprintf('%s.o', jobid));
  logerr       = fullfile(curPwd, sprintf('%s.e', jobid));
  fclose(fopen(logout, 'w'));
  fclose(fopen(logerr, 'w'));
end

% add the job to the persistent list, this is used for cleanup in case of Ctrl-C
qsublist('add', jobid, strtrim(result));

puttime = toc(stopwatch);

% remember the input arguments to speed up subsequent calls
previous_argin  = varargin;
