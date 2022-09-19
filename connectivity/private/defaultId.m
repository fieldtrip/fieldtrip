function id = defaultId(stack)

% DEFAULTID returns a string that can serve as warning or error identifier,
% for example 'FieldTip:ft_read_header:line345'.
%
% See also WARNING, ERROR, FT_NOTICE, FT_INFO, FT_DEBUG

% Copyright (C) 2017, Robert Oostenveld
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

if nargin==1
  % stack is already in input, and assumed to be properly pruned
else
  stack = dbstack('-completenames');
  
  % stack(1).name is defaultId
  % stack(2).name is ft_notification
  % stack(3).name is ft_error, ft_warning, ft_notice, ft_info or ft_debug
  
  % remove the functions that pertain to the notification system itself
  stack = stack(4:end);
  
  % remove the non-FieldTrip functions and scripts, these should not be part of the message identifier
  [v, p] = ft_version;
  keep   = startsWith({stack.file}, p);
  stack  = stack(keep);
end

switch numel(stack)
  case 0
    % it is called from the command line
    id = 'FieldTrip:commandline';
  case 1
    id = sprintf('FieldTrip:%s:line%d', stack(1).name, stack(1).line);
  case 2
    id = sprintf('FieldTrip:%s:%s:line%d', stack(2).name, stack(1).name, stack(1).line);
  otherwise
    % it is called from within a function
    name = sprintf('%s:', stack(end:-1:1).name); % this creates something like fun1:fun2:fun3:
    id   = sprintf('FieldTrip:%sline%d', name, stack(1).line);
end

% slashes occur when using nested functions, but are not allowed in the identifier
id(id=='/') = ':';
id(id=='\') = ':';
