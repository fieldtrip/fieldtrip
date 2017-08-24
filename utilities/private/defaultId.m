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
stack = stack(2:end); % this is this function itself

% remove the non-FieldTrip functions from the path, these should not be part of the default message identifier
keep = true(size(stack));
[v, p] = ft_version;
for i=1:numel(stack)
  keep(i) = strncmp(p, stack(i).file, length(p));
end
stack = stack(keep);


if ~isempty(stack)
  % it is called from within a function
  name = fliplr({stack.name});
  id = ['FieldTrip' sprintf(':%s', name{:}) ':line' num2str(stack(1).line)];
else
  % it is called from the command line
  id = 'FieldTrip:commandline';
end
