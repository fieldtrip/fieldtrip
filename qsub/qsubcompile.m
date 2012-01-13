function [fcomp] = qsubcompile(fname, batch)

% QSUBCOMPILE compiles your function into an fcomp that can easily be
% distributed on a cluster by QSUBCELLFUN.  Running a compiled version
% of your function does not take any additional MATLAB licenses. It does
% require that the appropriate MATLAB run-time environment is installed
% on your cluster.
%
% Use as
%   compiledfun = qsubcompile(fname)
%   argout      = qsubcellfun(compiledfun, argin, ...)
% or
%   compiledfun = qsubcompile(fname)
%   jobid       = qsubfexec(compiledfun, argin, ...)
%   argout      = qsubget(jobid)
%
% When executing a single batch of jobs using QSUBCELLFUN, you can also
% compile your function on the fly with the compile flag like this
%   argout      = qsubcellfun(fname, argin, ..., 'compile', 'yes')
% Using this syntax, the compiled function will be automatically cleaned
% up immediately after execution.
%
% A common problem for compilation is caused by the use of addpath in
% your startup.m file. Please change your startup.m file into
%   if ~isdeployed
%    % here goes the original content of your startup file
%    % ...
%   end
%
% See also QSUBCELLFUN, QSUBFEVAL, MCC, ISDEPLOYED

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

% Note that the function will be wrapped into qsubexec, which takes care of
% the loading and saving of the input and output arguments.

% It is possible for the user to have multiple MATLAB sessions running at
% the same time (on the same or different computers) and to start multiple
% instances of qsubcellfun. To ensure that those do not interfere with each
% others, each batch of jobs (i.e. instance of qsubcellfun) should get a
% unique identifier that is used in the filename of the temporary mat files.

if nargin<2 || isempty(batch)
  batch   = getbatch();              % this is a unique number
end
batchid = generatebatchid(batch);    % this is user_host_pid_batch

% some temporary filse are made during compilation, these flags determine
% whether they can be cleaned up afterwards
hasreadme = exist('./readme.txt', 'file');
hasmcclog = exist('./mccExcludedFiles.log', 'file');

if isa(fname, 'function_handle')
  % convert the function handle back into a string (e.g. @plus should be 'plus')
  fname = func2str(fname);
end

fprintf('compiling %s into %s\n', fname, batchid);
% try to compile into a stand-allone application
% ensure that cellfun is included, it might be needed for stacked jobs
mcc('-N', '-R', '-nodisplay', '-o', batchid, '-m', 'qsubexec', 'cellfun', fname);
fprintf('finished compiling\n');

if ~hasreadme
  delete('./readme.txt');
end

if ~hasmcclog
  delete('./mccExcludedFiles.log');
end

% reemmber all details
fcomp.fname       = fname;
fcomp.batch       = batch;
fcomp.batchid     = batchid;
fcomp.executable  = fullfile(pwd, sprintf('run_%s.sh', batchid));
