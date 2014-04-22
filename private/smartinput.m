function [newval, change] = smartinput(question, oldval)

% SMARTINPUT helper function for smart interactive input from the command line
%
% Use as
%   [newval, change] = smartinput(question, oldval)
%
% See also INPUT, PAUSE

% Copyright (C) 2006-2014, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

if ischar(oldval)
  newval = input(question, 's');
else
  newval = input(question);
end
if isempty(newval)
  newval = oldval;
  change = 0;
elseif isempty(oldval) && ~isempty(newval)
  change = 1;
elseif ischar(oldval) && strcmp(oldval, newval)
  change = 0;
elseif ~ischar(oldval) && all(size(oldval)==size(newval)) && all(oldval==newval)
  change = 0;
else
  change = 1;
end

