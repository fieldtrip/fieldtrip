function [s] = rmsubfield(s, f, v)

% RMSUBFIELD removes the contents of the specified field from a structure
% just like the standard Matlab RMFIELD function, except that you can also
% specify nested fields using a '.' in the fieldname. The nesting can be
% arbitrary deep.
%
% Use as
%   s = rmsubfield(s, 'fieldname')
% or as
%   s = rmsubfield(s, 'fieldname.subfieldname')
%
% See also SETFIELD, GETSUBFIELD, ISSUBFIELD

% Copyright (C) 2006-2013, Robert Oostenveld
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

if ~ischar(f)
  error('incorrect input argument for fieldname');
end

% remove the nested subfield using recursion
[t, f] = strtok(f, '.');
if any(f=='.')
  u = rmsubfield(getfield(s, t), f);
  s = setfield(s, t, u);
else
  s = rmfield(s, t);
end
