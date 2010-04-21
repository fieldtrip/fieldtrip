function [X,Y,Width,Height,Lbl] = read_lay(layoutname);

% READ_LAY reads an electrode or gradiometer layout file
% Layout files are used for topoplotting and multiplotting.
%
% Use as
%   [X, Y, Width, Height, Lbl] = read_lay(layoutname)

% Copyright (C) 2003, Ole Jensen
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

if ~exist(layoutname, 'file')
  error(sprintf('could not open layout file: %s', layoutname));
end

% discard the channel number
[chNum,X,Y,Width,Height,Lbl] = textread(layoutname,'%f %f %f %f %f %q');
