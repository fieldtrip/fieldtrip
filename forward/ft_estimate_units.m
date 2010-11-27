function unit = ft_estimate_units(size)

% FT_ESTIMATE_UNITS tries to determine the units of a geometrical object by
% looking at its size and by relating this to the size of the human
% brain.
%
% Use as
%   unit = ft_estimate_units(size)
%
% This function will return one of the following strings
%   'm'
%   'dm'
%   'cm'
%   'mm'

% Copyright (C) 2009, Robert Oostenveld
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
% $Id: ft_estimate_units.m 946 2010-04-21 17:51:16Z roboos $

% do some magic based on the size
unit = {'m', 'dm', 'cm', 'mm'};
indx = round(log10(size)+2-0.2);

if indx>length(unit)
  indx = length(unit);
  warning('assuming that the units are "%s"', unit{indx});
end

if indx<1
  indx = 1;
  warning('assuming that the units are "%s"', unit{indx});
end

unit = unit{indx};
