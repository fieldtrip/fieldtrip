function [r] = issubfield(s, f)

% ISSUBFIELD tests for the presence of a field in a structure just like the standard
% Matlab ISFIELD function, except that you can also specify nested fields
% using a '.' in the fieldname. The nesting can be arbitrary deep.
%
% Use as
%   f = issubfield(s, 'fieldname')
% or as
%   f = issubfield(s, 'fieldname.subfieldname')
%
% This function returns true if the field is present and false if the field
% is not present.
%
% See also ISFIELD, GETSUBFIELD, SETSUBFIELD

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

%try
%  getsubfield(s, f);    % if this works, then the subfield must be present
%  r = true;
%catch
%  r = false;                % apparently the subfield is not present
%end

t = textscan(f,'%s','delimiter','.');
t = t{1};
r = true;
for k = 1:numel(t)
  if isfield(s, t{k})
    s = s.(t{k});
  else
    r = false;
    return;
  end
end