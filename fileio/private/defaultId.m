function id = defaultId

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

stack = dbstack('-completenames');

% remove the functions that pertain to the notification system itself
keep = true(size(stack));
keep(strcmp({stack.name}, 'defaultId'))       = false; % this function itself
keep(strcmp({stack.name}, 'ft_notification')) = false; % this one is doing the work underneath ft_error/ft_warning/ft_notice/etc.
keep(strcmp({stack.name}, 'ft_error'))        = false;
keep(strcmp({stack.name}, 'ft_warning'))      = false;
keep(strcmp({stack.name}, 'ft_notice'))       = false;
keep(strcmp({stack.name}, 'ft_info'))         = false;
keep(strcmp({stack.name}, 'ft_debug'))        = false;
stack = stack(keep);

% remove the non-FieldTrip functions from the path, these should not be part of the default message identifier
keep = true(size(stack));
p = fileparts(mfilename('fullpath'));
% strip away '/utilities/private' where this function is located
p = p(1:end-18);
for i=1:numel(stack)
  keep(i) = strncmp(p, stack(i).file, length(p));
end
stack = stack(keep);

if ~isempty(stack)
  % it is called from within a function
  stack = flipud(stack);
  name  = {stack.name};
  id    = ['FieldTrip' sprintf(':%s', name{:}) ':line' num2str(stack(end).line)];
else
  % it is called from the command line
  id = 'FieldTrip:commandline';
end

% slashes occur when using nested functions, but are not allowed in the identifier
id(id=='/') = ':';
id(id=='\') = ':';

