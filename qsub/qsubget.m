function varargout = qsubget(jobid, varargin)

% QSUBGET get the output arguments after the remote job has been executed.
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
% Copyright (C) 2010, Robert Oostenveld
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

success = false;
while ~success && (timeout == 0 || toc(stopwatch)<timeout)
  
  % the code is largely shared with fieldtrip/peer/peerget.m
  % this section is the only part where it is different between peer and qsub
  
  curPwd = getcustompwd();
  inputfile    = fullfile(curPwd, sprintf('%s_input.mat', jobid));
  outputfile   = fullfile(curPwd, sprintf('%s_output.mat', jobid));
  logout       = fullfile(curPwd, sprintf('%s.o*', jobid)); % note the wildcard in the file name
  logerr       = fullfile(curPwd, sprintf('%s.e*', jobid)); % note the wildcard in the file name
  
  % the STDERR output log file should be empty, otherwise it indicates an error
  tmp = dir(logerr);
  if ~isempty(tmp) && tmp.bytes>0
    % show the error that was printed on STDERR
    type(fullfile(curPwd, tmp.name));
    error('the batch queue system returned an error for job %s, now aborting', jobid);
  end
  
  % the stdout and stderr log files are the last ones created
  % wait until they exist prior to reading the results
  if exist(outputfile, 'file') && isfile(logout) && isfile(logerr)
    % load the results from the output file
    tmp = load(outputfile);
    argout  = tmp.argout;
    options = tmp.optout;
    success = true;
    % clean up all temporary files
    % delete(inputfile); % this one has already been deleted in qsubexec immediately after loading it
    delete(outputfile);
    delete(logout);
    delete(logerr);
    % remove the job from the persistent list
    qsublist('del', jobid);
  elseif timeout == 0
    break; % only check once, no waiting here
  else
    % the job results have not arrived yet
    % wait a little bit and try again
    pause(sleep);
    continue
  end
  
end % while

if success
  
  % look at the optional arguments
  warn        = ft_getopt(options, 'lastwarn');
  err         = ft_getopt(options, 'lasterr');
  diarystring = ft_getopt(options, 'diary');
  
  fprintf('job %s returned, it required %s and %s\n', jobid, print_tim(ft_getopt(options, 'timused', nan)), print_mem(ft_getopt(options, 'memused', nan)));
  
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
  
  if strcmp(diary, 'error') && ~isempty(err)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('%% an error was detected inside MATLAB, the diary output of the remote execution follows \n');
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('%s', diarystring);
    closeline = true;
  elseif strcmp(diary, 'warning') && ~isempty(warn)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('%% a warning was detected inside MATLAB, the diary output of the remote execution follows\n');
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('%s', diarystring);
    closeline = true;
  elseif strcmp(diary, 'always')
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('%% the output of the remote execution follows\n');
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('%s', diarystring);
    closeline = true;
  else
    closeline = false;
  end
  if ~isempty(warn)
    warning(warn);
  end
  if ~isempty(err)
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
  if closeline
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function that detects a file, even with a wildcard in the filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = isfile(name)
tmp = dir(name);
status = length(tmp)==1 && ~tmp.isdir;

