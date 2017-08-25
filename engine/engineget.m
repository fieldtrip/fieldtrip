function varargout = engineget(jobid, varargin)

% ENGINEGET get the output arguments after the remote job has been executed.
%
% Use as
%   jobid  = enginefeval(fname, arg1, arg2, ...)
%   argout = engineget(jobid, ...)
%
% Optional arguments can be specified in key-value pairs and can include
%   StopOnError    = boolean (default = true)
%   timeout        = number, in seconds (default = 0, i.e. return immediately if output cannot be retrieved)
%   sleep          = number, in seconds (default = 0.01)
%   output         = string, 'varargout' or 'cell' (default = 'varargout')
%   diary          = string, can be 'always', 'warning', 'error' (default = 'error')
%
% See also ENGINEFEVAL, ENGINECELLFUN, ENGINEPOOL

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
%
% $Id$
% -----------------------------------------------------------------------

% the following are to speed up subsequent calls
persistent previous_varargin previous_output previous_diary previous_timeout previous_StopOnError

if isequal(previous_varargin, varargin)
  % prevent the ft_getopt function from being called, because it is slow
  % reuse the values from the previous call
  output      = previous_output;
  diary       = previous_diary;
  timeout     = previous_timeout;
  StopOnError = previous_StopOnError;
else
  % get the optional arguments
  output      = ft_getopt(varargin, 'output',      'varargout');
  diary       = ft_getopt(varargin, 'diary',       'error');
  StopOnError = ft_getopt(varargin, 'StopOnError', true);
  timeout     = ft_getopt(varargin, 'timeout',     5);
end

enghandle = enginepool('find', jobid);
if numel(enghandle)~=1
  error('FieldTrip:engine:jobNotFound', 'cannot locate the engine with this job');
end

success   = false;
stopwatch = tic;
while(toc(stopwatch)<timeout)
  if ~engine('isbusy', enghandle)
    argout  = engine('get', enghandle, 'argout');
    options = engine('get', enghandle, 'optout');
    engine('eval', enghandle, 'clear all', 1); % this should be a blocking call
    enginepool('release', enghandle);
    success = true;
    break
  else
    pause(0.1);
  end
end

if success
  % look at the optional arguments
  warn        = ft_getopt(options, 'lastwarn');
  err         = ft_getopt(options, 'lasterr');
  diarystring = ft_getopt(options, 'diary');
  
  fprintf('job %d returned, it required %s and %s\n', jobid, print_tim(ft_getopt(options, 'timused', nan)), print_mem(ft_getopt(options, 'memused', nan)));
  
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
  warning('FieldTrip:engine:jobNotAvailable', 'the job results are not yet available');
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
previous_output      = output;
previous_diary       = diary;
previous_timeout     = timeout;
previous_StopOnError = StopOnError;
