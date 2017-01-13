function [ws warned] = print_once(varargin)

% PRINT_ONCE will print a message for every unique point in the
% stacktrace only, e.g. in a for-loop a warning is thrown only once.
%
% Use as one of the following
%   print_once(string)
%   print_once(id, string)
% Alternatively, you can use ft_warning using a timeout
%   print_once(string, timeout)
%   print_once(id, string, timeout)
% where timeout should be inf if you don't want to see the warning ever
% again.
%
% Use as print_once('-clear') to clear old messages from the current
% stack
%
% Please note that you can NOT use it like this
%   print_once('the value is %d', 10)
% instead you should do
%   print_once(sprintf('the value is %d', 10))
%
% See also FT_WARNING

% Copyright (C) 2012, Robert Oostenveld
% Copyright (C) 2013, Robert Oostenveld, Jörn M. Horschig
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

if nargin < 1
  error('You need to specify at least a message');
end

warned = false;

if ~isfield(ft_default, 'printonce')
  ft_default.printonce.stopwatch  = [];
  ft_default.printonce.identifier = [];
end

% put the arguments we will pass to fprintf() in this cell array
warningArgs = {};

if nargin==3
  % calling syntax (id, msg, timeout)
  
  warningArgs = varargin(1:2);
  msg = warningArgs{2};
  timeout = varargin{3};
  fname = [warningArgs{1} '_' warningArgs{2}];
  
elseif nargin==2 && isnumeric(varargin{2})
  % calling syntax (msg, timeout)
  
  warningArgs = varargin(1);
  msg = warningArgs{1};
  timeout = varargin{2};
  fname = warningArgs{1};
  
elseif nargin==2 && ~isnumeric(varargin{2})
  % calling syntax (id, msg)
  
  warningArgs = varargin(1:2);
  msg = warningArgs{2};
  timeout = inf;
  fname = [warningArgs{1} '_' warningArgs{2}];
  
elseif nargin==1
  % calling syntax (msg)
  
  warningArgs = varargin(1);
  msg = warningArgs{1};
  timeout = inf; % default timeout in seconds
  fname = [warningArgs{1}];
  
end

if isempty(timeout)
  error('Timeout ill-specified');
end

if timeout ~= inf
  fname = decomma(fixname(fname)); % make a nice string that is allowed as structure fieldname
  if length(fname) > 63 % MATLAB max name
    fname = fname(1:63);
  end
  line = [];
else
  % here, we create the fieldname functionA.functionB.functionC... 
  [tmpfname ft_default.printonce.identifier line] = fieldnameFromStack(ft_default.printonce.identifier);
  if ~isempty(tmpfname)
    fname = tmpfname;
    clear tmpfname;
  end
end

if nargin==1 && ischar(varargin{1}) && strcmp('-clear', varargin{1})
  if strcmp(fname, '-clear') % reset all fields if called outside a function
    ft_default.printonce.identifier = [];
    ft_default.printonce.stopwatch  = [];
  else
    if issubfield(ft_default.printonce.identifier, fname)
      ft_default.printonce.identifier = rmsubfield(ft_default.printonce.identifier, fname);
    end
  end
  return;
end

% and add the line number to make this unique for the last function
fname = horzcat(fname, line);
  
if ~issubfield('ft_default.printonce.stopwatch', fname)
  ft_default.printonce.stopwatch = setsubfield(ft_default.printonce.stopwatch, fname, tic);
end

now = toc(getsubfield(ft_default.printonce.stopwatch, fname)); % measure time since first function call

if ~issubfield(ft_default.printonce.identifier, fname) || ...
    (issubfield(ft_default.printonce.identifier, fname) && now>getsubfield(ft_default.printonce.identifier, [fname '.timeout']))...
    || ~isequal(getsubfield(ft_default.printonce.identifier, [fname '.ws']), msg)

  % create or reset field
  ft_default.printonce.identifier = setsubfield(ft_default.printonce.identifier, fname, []);
    
  % message never given before or timed out, so proceed to print
  
  % escape characters
  escapedmsg = strrep(msg, '%', '%%');
  escapedmsg = strrep(escapedmsg, '\', '\\');
  fprintf(escapedmsg);
  
  % make sure ft_progress does not erase our message, in case it's running
  ft_default.progress.noerase = 1;
  
  ft_default.printonce.identifier = setsubfield(ft_default.printonce.identifier, [fname '.timeout'], now+timeout);
  ft_default.printonce.identifier = setsubfield(ft_default.printonce.identifier, [fname '.ws'], msg);
  warned = true;
else

  % the warning has been issued before, but has not timed out yet
  getsubfield(ft_default.printonce.identifier, [fname '.ws']);
  
end


end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function name = decomma(name)
name(name==',')=[];
end % function

function [fname ft_previous_warnings line] = fieldnameFromStack(ft_previous_warnings)
% stack(1) is this function, stack(2) is ft_warning
stack = dbstack('-completenames');
if size(stack) < 3
  fname = [];
  line = [];
  return;
end
i0 = 3;
% ignore ft_preamble
while strfind(stack(i0).name, 'ft_preamble')
  i0=i0+1;
end

fname = horzcat(fixname(stack(end).name));
if ~issubfield(ft_previous_warnings, fixname(stack(end).name))
  ft_previous_warnings.(fixname(stack(end).name)) = []; % iteratively build up structure fields
end
  

for i=numel(stack)-1:-1:(i0)
  % skip postamble scripts
  if strncmp(stack(i).name, 'ft_postamble', 12)
    break;
  end

  fname = horzcat(fname, '.', horzcat(fixname(stack(i).name))); % , stack(i).file
  if ~issubfield(ft_previous_warnings, fname) % iteratively build up structure fields
    setsubfield(ft_previous_warnings, fname, []);
  end
end

% line of last function call
line = ['.line', num2str(stack(i0).line)];
end

% function outcome = issubfield(strct, fname)
% substrindx = strfind(fname, '.');
% if numel(substrindx) > 0
%   % separate the last fieldname from all former
%   outcome = eval(['isfield(strct.' fname(1:substrindx(end)-1) ', ''' fname(substrindx(end)+1:end) ''')']);
% else
%   % there is only one fieldname
%   outcome = isfield(strct, fname);
% end
% end

% function strct = rmsubfield(strct, fname)
% substrindx = strfind(fname, '.');
% if numel(substrindx) > 0
%   % separate the last fieldname from all former
%   strct = eval(['rmfield(strct.' fname(1:substrindx(end)-1) ', ''' fname(substrindx(end)+1:end) ''')']);
% else
%   % there is only one fieldname
%   strct = rmfield(strct, fname);
% end
% end
