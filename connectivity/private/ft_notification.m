function [varargout] = ft_notification(varargin)

% FT_NOTIFICATION works mostly like the WARNING and ERROR commands in MATLAB and
% is called by FT_ERROR, FT_WARNING, FT_NOTICE, FT_INFO, FT_DEBUG. Note that you
% should not call this function directly.
%
% Some examples:
%  ft_info on
%  ft_info on msgId
%  ft_info off
%  ft_info off msgId
%  ft_info once
%  ft_info once msgId
%  ft_info on  backtrace
%  ft_info off backtrace
%  ft_info on  verbose
%  ft_info off verbose
%
%  ft_info query      % shows the status of all notifications
%  ft_info last       % shows the last notification
%  ft_info clear      % clears the status of all notifications
%  ft_info timeout 10 % sets the timeout (for 'once') to 10 seconds
%
% See also DEFAULTID

% Copyright (C) 2012-2017, Robert Oostenveld, J?rn M. Horschig
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

global ft_default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% determine who is calling

stack = dbstack('-completenames');
% st(1) is this function
% st(2) is the calling function
switch stack(2).name
  case 'ft_debug'
    level = 'debug';
  case 'ft_info'
    level = 'info';
  case 'ft_notice'
    level = 'notice';
  case 'ft_warning'
    level = 'warning';
  case 'ft_error'
    level = 'error';
  otherwise
    error('this function cannot be called from %s', stack(2).name);
end

% remove this function itself and the ft_xxx calling function
stack = stack(3:end);

% remove the non-FieldTrip functions from the path, these should not be part of the default message identifier
keep = true(size(stack));
p = fileparts(which('ft_defaults'));
for i=1:numel(stack)
  keep(i) = strncmp(p, stack(i).file, length(p));
end
stack = stack(keep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% handle the defaults

if isfield(ft_default, 'notification') && isfield(ft_default.notification, level)
  % this would be error, warning, notice, info, debug
  s = ft_default.notification.(level);
else
  s = [];
end

% start with an empty state structure
if isempty(s)
  s = struct('identifier', {}, 'state', {}, 'timestamp', {});
end

% set the default notification state
if ~ismember('all', {s.identifier})
  s = setstate(s, 'all', 'on');
end

% set the default backtrace state
defaultbacktrace = false;
if ~ismember('backtrace', {s.identifier})
  switch level
    case {'debug' 'info' 'notice'}
      s = setstate(s, 'backtrace', 'off');
    case 'warning'
      defaultbacktrace = true;
      t = warning('query', 'backtrace'); % get the default state
      s = setstate(s, 'backtrace', t.state);
    case 'error'
      s = setstate(s, 'backtrace', 'on');
  end % switch
end

% set the default verbose state
defaultverbose = false;
if ~ismember('verbose', {s.identifier})
  switch level
    case 'warning'
      defaultverbose = true;
      t = warning('query', 'verbose');% get the default state
      s = setstate(s, 'verbose', t.state);
    otherwise
      s = setstate(s, 'verbose', 'off');
  end
end

% set the default timeout
if ~ismember('timeout', {s.identifier})
  s = setstate(s, 'timeout', 60);
end

% set the last notification to empty
if ~ismember('last', {s.identifier})
  state.message    = '';
  state.identifier = '';
  state.stack      = struct('file', {}, 'name', {}, 'line', {});
  s = setstate(s, 'last', state);
end

if strcmp(level, 'warning')
  ws = warning;
  % warnings should be on in general
  warning on
  % the backtrace is handled by this function
  warning off backtrace
  % the verbose message is handled by this function
  warning off verbose
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the notification state according to the input

if numel(varargin)>0 && (isstruct(varargin{1}) || isempty(varargin{1}))
  ft_default.notification.(level) = varargin{1};
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% act according to the first input argument

if isempty(varargin)
  varargin{1} = 'query';
end

switch varargin{1}
  case 'on'
    if numel(varargin)>1
      % switch a specific item on
      msgId = varargin{2};
      s = setstate(s, msgId, 'on');
      if strcmp(msgId, 'backtrace')
        defaultbacktrace = false;
      end
      if strcmp(msgId, 'verbose')
        defaultverbose = false;
      end
      % return the message state of all
      varargout{1} = getreturnstate(s, msgId);
    else
      % switch all on
      s = setstate(s, 'all', 'on');
      % return the message state of all
      varargout{1} = getreturnstate(s);
    end
    
  case 'off'
    if numel(varargin)>1
      % switch a specific item off
      msgId = varargin{2};
      s = setstate(s, msgId, 'off');
      if strcmp(msgId, 'backtrace')
        defaultbacktrace = false;
      end
      if strcmp(msgId, 'verbose')
        defaultverbose = false;
      end
      % return the specific message state
      varargout{1} = getreturnstate(s, msgId);
    else
      % switch all off
      s = setstate(s, 'all', 'off');
      % return the message state of all
      varargout{1} = getreturnstate(s);
    end
    
  case 'once'
    if numel(varargin)>1
      % switch a specific item to once
      msgId = varargin{2};
      s = setstate(s, msgId, 'once');
      % return the specific message state
      varargout{1} = getreturnstate(s, msgId);
    else
      % switch all to once
      s = setstate(s, 'all', 'once');
      % return the message state of all
      varargout{1} = getreturnstate(s);
    end
    
  case 'timeout'
    % set the timeout, this is used for 'once'
    if ischar(varargin{2})
      s = setstate(s, 'timeout', str2double(varargin{2}));
    else
      s = setstate(s, 'timeout', varargin{2});
    end
    
  case {'last' '-last'}
    % return the last notification
    varargout{1} = getstate(s, 'last');
    
  case {'clear' '-clear'}
    % reset the notification system
    s = [];
    
  case {'query' '-query'}
    if numel(varargin)>1
      % select a specific item
      msgId = varargin{2};
      if ~ismember(msgId, {s.identifier})
        error('Unknown setting or incorrect message identifier ''%s''.', msgId);
      end
      msgState = getstate(s, msgId);
      if nargout
        varargout{1} = getreturnstate(s, msgId);
      elseif strcmp(msgId, 'verbose')
        if istrue(msgState)
          fprintf('%s output is verbose.\n', level);
        else
          fprintf('%s output is terse.\n', level);
        end
      elseif strcmp(msgId, 'backtrace')
        if istrue(msgState)
          fprintf('%s backtraces are enabled.\n', level);
        else
          fprintf('%s backtraces are disabled.\n', level);
        end
      else
        fprintf('The state of %s ''%s'' is ''%s''\n', level, msgId, msgState);
      end
    else
      % return all items
      r = getreturnstate(s);
      
      if nargout
        % return the state of all items
        varargout{1} = r;
      else
        % show the state of all items that are different from the default
        default = getstate(s, 'all');
        fprintf('The default %s state is ''%s''.', level, default);
        r = r(~strcmp({r.state}, default));
        % don't show these
        r(strcmp('verbose',   {r.identifier})) = [];
        r(strcmp('backtrace', {r.identifier})) = [];
        if ~isempty(r)
          fprintf(' Items not set to the default are\n\n');
        end
        for i=1:numel(r)
          fprintf('  %4s  %s\n', r(i).state, r(i).identifier);
        end
        fprintf('\n');
      end
    end
    
  otherwise
    
    if nargout
      varargout{1} = getreturnstate(s);
    end
    
    % first input might be msgId
    if numel(varargin)>1 && ~isempty(regexp(varargin{1}, '^\w*:', 'once'))
      msgId = varargin{1};
      varargin = varargin(2:end); % shift them all by one
    else
      msgId = defaultId;
    end
    
    % get the state for this notification, it will default to the 'all' state
    msgState = getstate(s, msgId);
    
    % errors are always to be printed
    if strcmp(level, 'error')
      msgState = 'on';
    end
    
    if strcmp(msgState, 'once')
      timeout = getstate(s, 'timeout');
      since   = elapsed(gettimestamp(s, msgId));
      if (since>timeout)
        % the timeout has passed, update the timestamp and print the message
        s = settimestamp(s, msgId, tic);
        msgState = 'on';
      else
        % the timeout has not yet passed, do not print the message
        msgState = 'off';
      end
    end
    
    % use an automatically generated default identifier
    if isempty(msgId)
      msgId = defaultId;
    end
    
    % store the last notification
    state.message    = strtrim(sprintf(varargin{:})); % remove the trailing newline
    state.identifier = msgId;
    state.stack      = stack;
    s = setstate(s, 'last', state);
    
    if strcmp(msgState, 'on')
      
      if strcmp(level, 'error')
        % update the global variable, we won't return here after the error
        ft_default.notification.(level) = s;
        % the remainder is fully handled by the ERROR function
        if ~isempty(msgId)
          error(msgId, varargin{:});
        else
          error(varargin{:});
        end
        
      elseif strcmp(level, 'warning')
        if ~isempty(msgId)
          warning(msgId, varargin{:});
        else
          warning(varargin{:});
        end
        
      else
        % ensure there is a line end
        if isempty(regexp(varargin{1}, '\\n$', 'once')) % note the double \\
          varargin{1} = [varargin{1} '\n'];
        end
        fprintf(varargin{:});
        
      end % if level=error, warning or otherwise
      
      % decide whether the stack trace should be shown
      if istrue(getstate(s, 'backtrace'))
        for i=1:numel(stack)
          % show the deepest and lowest-level function first
          [p, f, x] = fileparts(stack(i).file);
          if isequal(f, stack(i).name)
            funname = stack(i).name;
          else
            funname = sprintf('%s>%s', f, stack(i).name);
          end
          filename = stack(i).file;
          if isempty(fileparts(filename))
            % it requires the full path
            filename = fullfile(pwd, filename);
          end
          if ft_platform_supports('html')
            fprintf(' In <a href = "matlab: opentoline(''%s'',%d,1)">%s at line %d</a>\n', filename, stack(i).line, funname, stack(i).line);
          else
            fprintf(' In ''%s'' at line %d\n', filename, stack(i).line);
          end
        end
        fprintf('\n');
      end
      
      % decide whether a verboe message should be shown
      if istrue(getstate(s, 'verbose'))
        if ~isempty(msgId)
          fprintf('Type "ft_%s off %s" to suppress this message.\n', level, msgId)
        else
          fprintf('Type "ft_%s off" to suppress this message.\n', level)
        end
      end
      
    end % if msgState is on
    
end % switch varargin{1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update the global variable
if defaultbacktrace && ~isempty(s)
  % do not keep the default, it will be determined again on the next call
  s = s(~strcmp({s.identifier}, 'backtrace'));
end
if defaultverbose && ~isempty(s)
  % do not keep the default, it will be determined again on the next call
  s = s(~strcmp({s.identifier}, 'verbose'));
end
ft_default.notification.(level) = s;

if strcmp(level, 'warning')
  % return to the original warning state
  warning(ws);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = getreturnstate(s, msgId)
if nargin<2
  r = s;
  % don't return these
  r(strcmp('timeout',   {r.identifier})) = [];
  r(strcmp('last',      {r.identifier})) = [];
  % don't return the timestamps
  r = rmfield(r, 'timestamp');
else
  msgState = getstate(s, msgId);
  r = struct('identifier', msgId, 'state', msgState);
end

function state = getstate(s, msgId)
identifier = {s.identifier};
sel = find(strcmp(identifier, msgId));
if numel(sel)==1
  state = s(sel).state;
  if isempty(state)
    state = getstate(s, 'all');
  end
else
  state = getstate(s, 'all');
end

function timestamp = gettimestamp(s, msgId)
identifier = {s.identifier};
sel = find(strcmp(identifier, msgId));
if numel(sel)==1
  timestamp = s(sel).timestamp;
else
  timestamp = nan;
end

function s = setstate(s, msgId, state)
if isempty(msgId), return; end % this happens from the command line
identifier = {s.identifier};
sel = find(strcmp(identifier, msgId));
if numel(sel)==1
  s(sel).state = state;
else
  s(end+1).identifier = msgId;
  s(end  ).state      = state;
  s(end  ).timestamp  = nan;
end

function s = settimestamp(s, msgId, timestamp)
if isempty(msgId), return; end % this happens from the command line
identifier = {s.identifier};
sel = find(strcmp(identifier, msgId));
if numel(sel)==1
  s(sel).timestamp = timestamp;
else
  s(end+1).identifier = msgId;
  s(end  ).state      = [];
  s(end  ).timestamp  = timestamp;
end

function t = elapsed(start)
if isempty(start) || isnan(start)
  t = inf;
else
  t = toc(start);
end
