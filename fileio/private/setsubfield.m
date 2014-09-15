function [s] = setsubfield(s, f, v)

% SETSUBFIELD sets the contents of the specified field to a specified value
% just like the standard Matlab SETFIELD function, except that you can also
% specify nested fields using a '.' in the fieldname. The nesting can be
% arbitrary deep.
%
% Use as
%   s = setsubfield(s, 'fieldname', value)
% or as
%   s = setsubfield(s, 'fieldname.subfieldname', value)
%
% where nested is a logical, false denoting that setsubfield will create
% s.subfieldname instead of s.fieldname.subfieldname
%
% See also SETFIELD, GETSUBFIELD, ISSUBFIELD

% Copyright (C) 2005-2013, Robert Oostenveld
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

% t = {};
% while (1)
%   [t{end+1}, f] = strtok(f, '.');
%   if isempty(f)
%     break
%   end
%   
% end
%
%s = setfield(s, t{:}, v);

t = textscan(f,'%s','delimiter','.');
s = setfield(s, t{1}{:}, v);
