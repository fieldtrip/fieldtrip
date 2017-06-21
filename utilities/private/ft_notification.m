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
%  ft_info backtrace on
%  ft_info backtrace off
%
%  ft_info query      % shows the status of all notifications
%  ft_info last       % shows the last notification
%  ft_info clear      % clears the status of all notifications
%  ft_info timeout 10 % sets the timeout (for 'once') to 10 seconds

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

% remove this function itself and the calling function
stack = stack(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the notification state according to the input

if numel(varargin)>0 && isstruct(varargin{1})
  ft_default.notification.(level) = varargin{1};
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% handle the defaults

if isfield(ft_default, 'notification') && isfield(ft_default.notification, level)
  % this would be error, warning, notice, info, debug
  s = ft_default.notification.(level);
else
  s = [];
end

% set the default notification state
if isempty(s) || ~ismember('all', {s.identifier})
  s(end+1).identifier = 'all';
  s(end  ).state      = 'on';
  s(end  ).timestamp  = nan;  % this ensures that the field exists
end

% set the default backtrace state
if isempty(s) || ~ismember('backtrace', {s.identifier})
  s(end+1).identifier = 'backtrace';
  s(end  ).state      = 'on';
  s(end  ).timestamp  = nan;  % this ensures that the field exists
end

% set the default timeout
if isempty(s) || ~ismember('timeout', {s.identifier})
  s(end+1).identifier = 'timeout';
  s(end  ).state      = 60;   % default is 60 seconds
  s(end  ).timestamp  = nan;  % this ensures that the field exists
end

% set the last notification to empty
if isempty(s) || ~ismember('last', {s.identifier})
  s(end+1).identifier = 'last';
  s(end  ).state.message    = '';
  s(end  ).state.identifier = '';
  s(end  ).state.stack      = struct('file', {}, 'name', {}, 'line', {});
  s(end  ).timestamp  = nan;  % this ensures that the field exists
end

if isempty(varargin)
  varargin{1} = 'query';
end

if ~isempty(stack)
  % it is called from within a function
  name = {stack.name};
  defaultId = ['FieldTrip' sprintf(':%s', name{:}) ':' num2str(stack(1).line)];
else
  % it is called from the command line
  defaultId = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% act according to the first input argument

switch varargin{1}
  case 'on'
    if numel(varargin)>1
      % switch a specific item on
      msgId = varargin{2};
      s = setstate(s, msgId, 'on');
    else
      % switch all on
      s = setstate(s, 'all', 'on');
    end
    
  case 'off'
    if numel(varargin)>1
      % switch a specific item off
      msgId = varargin{2};
      s = setstate(s, msgId, 'off');
    else
      % switch all off
      s = setstate(s, 'all', 'off');
    end
    
  case 'once'
    if numel(varargin)>1
      % switch a specific item to once
      msgId = varargin{2};
      s = setstate(s, msgId, 'once');
    else
      % switch all to once
      s = setstate(s, 'all', 'once');
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
      msgState = getstate(s, msgId);
      if nargout
        r = struct('identifier', msgId, 'state', msgState);
        varargout{1} = r;
      else
        fprintf('The state of %s ''%s'' is ''%s''\n', level, msgId, msgState);
      end
    else
      % return all items
      r = s;
      % don't return the backtrace, timeout and last
      r(strcmp('backtrace', {r.identifier})) = [];
      r(strcmp('timeout', {r.identifier})) = [];
      r(strcmp('last', {r.identifier})) = [];
      
      if nargout
        % do not return the timestamp field
        varargout{1} = rmfield(r, 'timestamp');
      else
        % show the state of all items that are different from the default
        default = getstate(s, 'all');
        fprintf('The default %s state is ''%s''.', level, default);
        r = r(~strcmp({r.state}, default));
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
    % first input might be msgId
    if any(varargin{1}==':') && numel(varargin)>1
      msgId = varargin{1};
      varargin = varargin(2:end); % shift them all by one
    else
      msgId = [];
    end
    
    % get the state for this notification, it will default to the 'all' state
    msgState = getstate(s, msgId);
    
    % errors are always to be printed
    if strcmp(level, 'error')
      msgState = 'on';
    end
    
    % ensure there is a line end
    if isempty(regexp(varargin{1}, '\\n$', 'once')) % note the double \\
      varargin{1} = [varargin{1} '\n'];
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
    state.message    = sprintf(varargin{:});
    state.message    = state.message(1:end-1); % remove the trailing newline
    state.identifier = msgId;
    state.stack      = stack;
    s = setstate(s, 'last', state);
    
    if strcmp(msgState, 'on')
      
      if strcmp(level, 'error')
        % update the global variable
        ft_default.notification.(level) = s;
        % the remainder is fully handled by the ERROR function
        error(msgId, varargin{:});
        
      elseif strcmp(level, 'warning')
        % the backtrace is handled by the WARNING function
        if istrue(getstate(s, 'backtrace'))
          warning('backtrace', 'on');
        else
          warning('backtrace', 'off');
        end
        warning(msgId, varargin{:});
        
      else
        fprintf(varargin{:});
        
        % decide whether the stack should be shown
        if istrue(getstate(s, 'backtrace'))
          for i=numel(stack):-1:1
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
            fprintf(' In <a href = "matlab: opentoline(''%s'',%d,1)">%s at line %d</a>\n', filename, stack(i).line, funname, stack(i).line);
          end
          fprintf('\n');
        end
        
      end % if level=error, warning or otherwise
      
    else
      % don't print it
    end % if msgState is on
    
end % switch varargin{1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update the global variable
ft_default.notification.(level) = s;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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