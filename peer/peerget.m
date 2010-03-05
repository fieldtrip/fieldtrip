function varargout = peerget(jobid, varargin)

% PEERGET get the output arguments after the remote job has been executed.
%
% Use as
%   argout = peerget(jobid, ...)
%
% Optional arguments can be specified in key-value pairs and can include
%   timeout  = number, in seconds (default = 1)
%   sleep    = number, in seconds (default = 0.01)
%   output   = string, 'varargout' or 'cell' (default = 'varargout')
%   diary    = string, can be 'always', 'warning', 'error' (default = 'error')
%
% See also PEERFEVAL, PEERCELLFUN

% get the optional arguments
timeout = keyval('timeout', varargin); if isempty(timeout), timeout=1;          end
sleep   = keyval('sleep',   varargin); if isempty(sleep),   sleep=0.01;         end
output  = keyval('output',  varargin); if isempty(output),  output='varargout'; end
diary   = keyval('diary',   varargin); if isempty(diary),   diary='error';      end

% keep track of the time
stopwatch = tic;

success = false;
while ~success && toc(stopwatch)<timeout

  joblist = peer('joblist');
  sel = find([joblist.jobid]==jobid);

  if ~isempty(sel)
    [argout, options] = peer('get', jobid);
    peer('clear', jobid);
    success = true;
  else
    % the job results have not arrived yet
    % wait a little bit and try again
    pause(sleep);
    continue
  end

end % while 

if success

  % look at the optional arguments
  elapsed     = keyval('elapsed',  options);
  warn        = keyval('lastwarn', options);
  err         = keyval('lasterr',  options);
  diarystring = keyval('diary',    options);

  if strcmp(diary, 'error') && ~isempty(err)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('%% an error was detected, the diary output of the remote execution follows \n');
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('%s', diarystring);
    closeline = true;
  elseif strcmp(diary, 'warning') && ~isempty(warn)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('%% a warning was detected, the diary output of the remote execution follows\n');
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
    if ischar(err)
      % it only contains the description
      error(err);
    else
      % it contains the full details
      ws = warning('off', 'MATLAB:structOnObject');
      rethrow(struct(err));
      warning(ws);
    end
  end
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
  warning('the job results are not yet available');
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


