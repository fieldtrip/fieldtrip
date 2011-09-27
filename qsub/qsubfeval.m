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

% these are used to speed up the processing of multiple function calls with
% the same input arguments (e.g. from peercellfun)
persistent previous_argin

% keep track of the time
stopwatch = tic;

% convert the input arguments into something that strmatch can work with
strargin = varargin;
strargin(~cellfun(@ischar, strargin)) = {''};

% locate the begin of the optional key-value arguments
optbeg = false(size(strargin));
optbeg = optbeg | strcmp('memreq',  strargin);
optbeg = optbeg | strcmp('timreq',  strargin);
optbeg = optbeg | strcmp('diary',   strargin);
optbeg = optbeg | strcmp('batch',    strargin);
optbeg = optbeg | strcmp('timoverhead', strargin);
optbeg = optbeg | strcmp('memoverhead', strargin);
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

% skip the optional key-value arguments
if ~isempty(optbeg)
  varargin = varargin(1:(optbeg-1));
end

if isa(varargin{1}, 'function_handle')
  % convert the function handle back into a string (e.g. @plus should be 'plus')
  varargin{1} = func2str(varargin{1});
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
  matlabcmd = 'matlab78 -singleCompThread';
elseif matlabversion(7.9) % 2009b
  matlabcmd = 'matlab79 -singleCompThread';
elseif matlabversion('2010a')
  matlabcmd = 'matlab2010a -singleCompThread';
elseif matlabversion('2010b')
  matlabcmd = 'matlab2010b -singleCompThread';
elseif matlabversion('2011a')
  matlabcmd = 'matlab2011a -singleCompThread';
elseif matlabversion('2011b')
  matlabcmd = 'matlab2011b -singleCompThread';
elseif matlabversion('2012a')
  matlabcmd = 'matlab2012a -singleCompThread';
elseif matlabversion('2012b')
  matlabcmd = 'matlab2012b -singleCompThread';
else
  % use whatever is available as default
  matlabcmd = 'matlab';
end

% create the matlab script commands (one entry per line)
matlabScript = [...
  'restoredefaultpath;',...
  sprintf('addpath(''%s'');', fileparts(mfilename('fullpath'))),...
  sprintf('qsubexec(''%s'');', jobid),...
  sprintf('exit')];

% create the shell commands to execute matlab (one entry per line)
shellCmd = [...
  sprintf('cd \\"%s\\"\n', curPwd),...
  sprintf('%s -nosplash -nodisplay -r \\"%s\\"\n', matlabcmd, matlabScript)];

% set the job requirements according to the users specification
requirements = '';
if ~isempty(timreq)
  requirements = [requirements sprintf('-l walltime=%d ', timreq+timoverhead)];
end
if ~isempty(memreq)
  % don't know the difference
  requirements = [requirements sprintf('-l mem=%.0f ',   memreq+memoverhead)];
%   requirements = [requirements sprintf('-l vmem=%.0f ',  memreq+memoverhead)];
%   requirements = [requirements sprintf('-l pmem=%.0f ',  memreq+memoverhead)];
%   requirements = [requirements sprintf('-l pvmem=%.0f ', memreq+memoverhead)];
end

% In the command below both stderr and stout are redirected to /dev/null, so any
% output information will not be available for inspection. However, any
% matlab errors will be reported back by fexec.
% cmdline = ['qsub -e /dev/null -o /dev/null -N ' jobid ' ' requirements shellscript];

% qsubfget will check the stderr output log file for errors, e.g. MATLAB crashes 
cmdline = sprintf('echo "%s" | qsub -N %s %s -o %s -e %s', ...
  shellCmd, jobid, requirements, curPwd, curPwd);

fprintf('submitting job %s...', jobid); 
[status,result] = system(cmdline);
fprintf(' qstat job id %s\n', strtrim(result));

% add the job to the persistent list, this is used for cleanup in case of ctrl-C
qsublist('add', jobid, strtrim(result));

puttime = toc(stopwatch);

% remember the input arguments to speed up subsequent calls
previous_argin  = varargin;

