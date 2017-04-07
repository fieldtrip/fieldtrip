function test_xunit_stacktest

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_defaults ft_warning

ft_defaults

stack = dbstack('-completenames');
if size(stack) < 1
  fname = [];
  line = [];
  return;
end
i0 = 1;

fname = horzcat(stack(end).name);

fprintf('The fieldname would be tmpstrct.%s\n', fname);

tmpstrct = [];
setsubfield(tmpstrct, fname, [])
tmpstrct2 = [];
eval(['tmpstrct2.' fname ' =[]'])

for i=1:10
    ft_warning('Tthis is a test');
end
ft_warning('Tthis is another test');

ft_warning('-clear');

fprintf('Completed!\n');

end


function [ws warned] = ft_warning(varargin)
%
% WARNING_ONCE will throw a warning for every unique point in the
% stacktrace only, e.g. in a for-loop a warning is thrown only once.
%
% Use as one of the following
%   ft_warning(string)
%   ft_warning(id, string)
% Alternatively, you can use ft_warning using a timeout
%   ft_warning(string, timeout)
%   ft_warning(id, string, timeout)
% where timeout should be inf if you don't want to see the warning ever
% again.
%
% Use as ft_warning('-clear') to clear old warnings from the current
% stack
%
% It can be used instead of the MATLAB built-in function WARNING, thus as
%   s = ft_warning(...)
% or as
%   ft_warning(s)
% where s is a structure with fields 'identifier' and 'state', storing the
% state information. In other words, ft_warning accepts as an input the
% same structure it returns as an output. This returns or restores the
% states of warnings to their previous values.
%
% It can also be used as
%    [s w] = ft_warning(...)
% where w is a boolean that indicates whether a warning as been thrown or not.
%
% Please note that you can NOT use it like this
%   ft_warning('the value is %d', 10)
% instead you should do
%   ft_warning(sprintf('the value is %d', 10))
%

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

fprintf('Dummy ft_warning output - this is the new function\n');

if nargin < 1
  error('You need to specify at least a warning message');
end

warned = false;
if isstruct(varargin{1})
  warning(varargin{1});
  return;
end

if ~isfield(ft_default, 'warning')
  ft_default.warning.stopwatch  = [];
  ft_default.warning.identifier = [];
end

% put the arguments we will pass to warning() in this cell array
warningArgs = {};

if nargin==3
  % calling syntax (id, msg, timeout)
  
  warningArgs = varargin(1:2);
  timeout = varargin{3};
  fname = [warningArgs{1} '_' warningArgs{2}];
  
elseif nargin==2 && isnumeric(varargin{2})
  % calling syntax (msg, timeout)
  
  warningArgs = varargin(1);
  timeout = varargin{2};
  fname = warningArgs{1};
  
elseif nargin==2 && ~isnumeric(varargin{2})
  % calling syntax (id, msg)
  
  warningArgs = varargin(1:2);
  timeout = inf;
  fname = [warningArgs{1} '_' warningArgs{2}];
  
elseif nargin==1
  % calling syntax (msg)
  
  warningArgs = varargin(1);
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
  [tmpfname ft_default.warning.identifier line] = fieldnameFromStack(ft_default.warning.identifier);
  if ~isempty(tmpfname)
    fname = tmpfname;
    clear tmpfname;
  end
end

if nargin==1 && ischar(varargin{1}) && strcmp('-clear', varargin{1})
  if strcmp(fname, '-clear') % reset all fields if called outside a function
    ft_default.warning.identifier = [];
    ft_default.warning.stopwatch  = [];
  else
    ft_default.warning.identifier = rmsubfield(ft_default.warning.identifier, fname);
  end
  return;
end

% and add the line number to make this unique for the last function
fname = horzcat(fname, line);
  
if ~issubfield('ft_default.warning.stopwatch', fname)
  ft_default.warning.stopwatch = setsubfield(ft_default.warning.stopwatch, fname, tic);
end

now = toc(getsubfield(ft_default.warning.stopwatch, fname)); % measure time since first function call

if ~issubfield(ft_default.warning.identifier, fname) || ...
    (issubfield(ft_default.warning.identifier, fname) && now>getsubfield(ft_default.warning.identifier, [fname '.timeout']))

  % create or reset field
  ft_default.warning.identifier = setsubfield(ft_default.warning.identifier, fname, []);
    
  % warning never given before or timed out
  ws = warning(warningArgs{:});
  ft_default.warning.identifier = setsubfield(ft_default.warning.identifier, [fname '.timeout'], now+timeout);
  ft_default.warning.identifier = setsubfield(ft_default.warning.identifier, [fname '.ws'], ws);
  warned = true;
else

  % the warning has been issued before, but has not timed out yet
  ws = getsubfield(ft_default.warning.identifier, [fname '.ws']);
  
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

fname = horzcat(stack(end).name);
if ~issubfield(ft_previous_warnings, stack(end).name)
  ft_previous_warnings.(stack(end).name) = []; % iteratively build up structure fields
end
  

for i=numel(stack)-1:-1:(i0)
  fname = horzcat(fname, '.', horzcat(stack(i).name)); % , stack(i).file
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
