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
optbeg = optbeg | strcmp('timeout', strargin);
optbeg = optbeg | strcmp('sleep',   strargin);
optbeg = optbeg | strcmp('memreq',  strargin);
optbeg = optbeg | strcmp('cpureq',  strargin);
optbeg = optbeg | strcmp('timreq',  strargin);
optbeg = optbeg | strcmp('hostid',  strargin);
optbeg = optbeg | strcmp('diary',   strargin);
optbeg = find(optbeg);
optarg = varargin(optbeg:end);

% get the optional input arguments
timeout = ft_getopt(optarg, 'timeout', inf);
sleep   = ft_getopt(optarg, 'sleep',   0.05);
memreq  = ft_getopt(optarg, 'memreq',  0);
cpureq  = ft_getopt(optarg, 'cpureq',  0);
timreq  = ft_getopt(optarg, 'timreq',  0);
hostid  = ft_getopt(optarg, 'hostid',  []);
diary   = ft_getopt(optarg, 'diary',   []);

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

jobid = round(rand(1)*1e8);

% each job should have a different random number sequence
randomseed = rand(1)*double(intmax);

% pass some options that influence the remote execution
options = {'pwd', getcustompwd, 'path', getcustompath, 'global', getglobal, 'diary', diary, 'memreq', memreq, 'cpureq', cpureq, 'timreq', timreq, 'randomseed', randomseed};

p = getenv('HOME');
inputfile    = fullfile(p, sprintf('job_%08d_input.mat', jobid));
shellscript  = fullfile(p, sprintf('job_%08d.sh', jobid));
matlabscript = fullfile(p, sprintf('job_%08d.m', jobid));

% rename and save the variables
argin = varargin;
optin = options;
save(inputfile, 'argin', 'optin');

% create the shell script
fid = fopen(shellscript, 'wt');
fprintf(fid, '#!/bin/sh\n');
fprintf(fid, 'cd "%s"\n', p);
fprintf(fid, 'matlab2010b -nosplash -nodisplay -r job_%08d\n', jobid);
fclose(fid);

% create the matlab script
fid = fopen(matlabscript, 'wt');
fprintf(fid, 'restoredefaultpath\n');
fprintf(fid, 'addpath %s\n', fileparts(mfilename('fullpath')));
fprintf(fid, 'qsubexec(%d)\n', jobid);
fprintf(fid, 'exit\n');
fclose(fid);

cmdline = sprintf('qsub %s', shellscript);
% fprintf('submitting job %08d\n', jobid); 
system(cmdline);
puttime = toc(stopwatch);

% remember the input arguments to speed up subsequent calls
previous_argin  = varargin;

