function varargout = qsubget(jobid, varargin)

% QSUBGET get the output arguments after the remote job has been executed
% on the Torque, SGE, PBS or SLURM batch queue system.
%
% Use as
%   jobid  = qsubfeval(fname, arg1, arg2, ...)
%   argout = qsubget(jobid, ...)
%
% Optional arguments can be specified in key-value pairs and can include
%   StopOnError    = boolean (default = true)
%   timeout        = number, in seconds (default = 0; return immediately if output cannot be gotten)
%   sleep          = number, in seconds (default = 0.01)
%   output         = string, 'varargout' or 'cell' (default = 'varargout')
%   diary          = string, can be 'always', 'warning', 'error' (default = 'error')
%
% See also QSUBFEVAL, QSUBCELLFUN

% -----------------------------------------------------------------------
% Copyright (C) 2010-2016, Robert Oostenveld
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

% the following are to speed up subsequent calls
persistent previous_varargin previous_timeout previous_sleep previous_output previous_diary previous_StopOnError

if isequal(previous_varargin, varargin)
  % prevent the ft_getopt function from being called, because it is slow
  % reuse the values from the previous call
  timeout     = previous_timeout;
  sleep       = previous_sleep;
  output      = previous_output;
  diary       = previous_diary;
  StopOnError = previous_StopOnError;
else
  % get the optional arguments
  timeout     = ft_getopt(varargin, 'timeout',     0);
  sleep       = ft_getopt(varargin, 'sleep',       0.010);
  output      = ft_getopt(varargin, 'output',      'varargout');
  diary       = ft_getopt(varargin, 'diary',       'error');
  StopOnError = ft_getopt(varargin, 'StopOnError', true);
end

% keep track of the time
stopwatch = tic;
completed = false;

while ~completed && (timeout == 0 || toc(stopwatch)<timeout)
  % the code is largely shared with fieldtrip/peer/peerget.m
  % this section is the only part where it is different between peer and qsub
  
  if qsublist('completed', jobid)
    % the job has completed and the output file exist
    
    curPwd       = getcustompwd();
    outputfile   = fullfile(curPwd, sprintf('%s_output.mat', jobid));
    logout       = fullfile(curPwd, sprintf('%s.o*', jobid)); % note the wildcard in the file name
    logerr       = fullfile(curPwd, sprintf('%s.e*', jobid)); % note the wildcard in the file name
    
    if exist(outputfile, 'file')
      tmp = load(outputfile);
      argout  = tmp.argout;
      options = tmp.optout;
    else
      argout  = {};
      options = {};
    end
    
    % the stderr log file should be empty, otherwise it potentially indicates an error
    direrr = dir(logerr);
    if ~exist(outputfile, 'file') && ~isempty(direrr) && direrr.bytes>0
      fid = fopen(fullfile(curPwd, direrr.name), 'r');
      stderr = fread(fid, [1 inf], 'char');
      fclose(fid);
    else
      stderr = '';
    end
    
    % the stdout log file can be used instead of the MATLAB diary
    dirout = dir(logout);
    if ~exist(outputfile, 'file') && ~isempty(dirout) && dirout.bytes>0
      fid = fopen(fullfile(curPwd, dirout.name), 'r');
      stdout = fread(fid, [1 inf], 'char');
      fclose(fid);
    else
      stdout = '';
    end
    
    % clean up all temporary files
    % delete(inputfile); % this one has already been deleted in qsubexec immediately after loading it
    if exist(outputfile, 'file'), delete(outputfile); end
    if ~isempty(dir(logout)), delete(logout); end % note the wildcard in the file name
    if ~isempty(dir(logerr)), delete(logerr); end % note the wildcard in the file name
    
    % remove the job from the persistent list
    qsublist('del', jobid);
    
    % indicate that the job has completed
    completed = true;
    
  elseif timeout == 0
    % only check once, no waiting here
    break;
    
  else
    % the job results have not arrived yet
    % wait a little bit and try again
    pausejava(sleep);
    continue
  end
  
end % while not completed

if completed
  
  if isempty(options)
    % the job did not complete normally, it might be aborted by the batch queuing system
    if ~isempty(stderr)
      % the stderr output can contain details
      err = sprintf('Error in job execution, the output on stderr is\n%s', stderr);
    else
      err = sprintf('Error in job execution');
    end
    % use the stderr and stdout output rather than the MATLAB results
    options = {'lasterr', err, 'diary', stdout};
  elseif isempty(ft_getopt(options, 'lasterr'))
    fprintf('job %s returned, it required %s and %s on %s\n', jobid, print_tim(ft_getopt(options, 'timused', nan)), print_mem(ft_getopt(options, 'memused', nan)), ft_getopt(options, 'hostname', 'unknown'));
  end
  
  % look at the optional arguments
  warn        = ft_getopt(options, 'lastwarn');
  err         = ft_getopt(options, 'lasterr');
  diarystring = ft_getopt(options, 'diary');
  
  % if there is an error, it needs to be represented as a message string
  % and optionally also as a strucure for rethrowing
  if ~isempty(err)
    if ischar(err)
      errmsg = err;
    elseif isstruct(err)
      errmsg = err.message;
    else
      errmsg = err.message;
      % convert the MEexception object into a structure to allow a rethrow further down in the code
      ws = warning('off', 'MATLAB:structOnObject');
      err = struct(err);
      warning(ws);
    end
  end
  
  % the information that is printed on screen starts and ends with a separator line
  separatorline = false;
  
  if ~isempty(diarystring)
    if strcmp(diary, 'error') && ~isempty(err)
      fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
      fprintf('%% an error was detected, the output of the remote execution follows \n');
      fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
      fprintf('%s', diarystring);
      separatorline = true;
    elseif strcmp(diary, 'warning') && ~isempty(warn)
      fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
      fprintf('%% a warning was detected, the output of the remote execution follows\n');
      fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
      fprintf('%s', diarystring);
      separatorline = true;
    elseif strcmp(diary, 'always')
      fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
      fprintf('%% no problem was detected, the output of the remote execution follows\n');
      fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
      fprintf('%s', diarystring);
      separatorline = true;
    end
  end
  
  if ~isempty(warn)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    warning(warn);
  end
  
  if ~isempty(err)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    if StopOnError
      if ischar(err)
        error(err);
      else
        rethrow(err);
      end
    else
      warning('error during remote execution: %s', errmsg);
    end
  end % ~isempty(err)
  
  if separatorline
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  end
  
  switch output
    case 'varargout'
      % return the output arguments, the options cannot be returned
      varargout = argout;
    case 'cell'
      % return the output arguments and the options
      varargout{1} = argout;
      varargout{2} = options;
    otherwise
      error('invalid output option');
  end
  
else
  warning('FieldTrip:qsub:jobNotAvailable', 'the job results are not yet available');
  switch output
    case 'varargout'
      % return empty output arguments
      varargout = cell(1, nargout);
    case 'cell'
      % return the output arguments and the options as empty cells
      varargout{1} = {};
      varargout{2} = {};
    otherwise
      error('invalid output option');
  end
end

% remember the input arguments to speed up subsequent calls
previous_varargin    = varargin;
previous_timeout     = timeout;
previous_sleep       = sleep;
previous_output      = output;
previous_diary       = diary;
previous_StopOnError = StopOnError;
