function y = istrue(x)

% ISTRUE converts an input argument like "yes/no", "true/false" or "on/off" into a
% boolean. If the input is boolean, then it will remain like that.

% Copyright (C) 2009-2012, Robert Oostenveld
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

true_list  = {'yes' 'true' 'on' 'y' };
false_list = {'no' 'false' 'off' 'n' 'none'};

if ischar(x)
  % convert string to boolean value
  if any(strcmpi(x, true_list))
    y = true;
  elseif any(strcmpi(x, false_list))
    y = false;
  else
    error('cannot determine whether "%s" should be interpreted as true or false', x);
  end
else
  % convert numerical value to boolean
  y = logical(x);
end

