function [h] = crosshair(pos, varargin)

% CROSSHAIR adds a crosshair at position (x,y) to the current plot
% additional options are passed to the builtin line function
% the handles of the lines are returned
% 
% h = crosshair([x,y])
% 
% see also: LINE, TEXT 

% Copyright (C) 2003, Robert Oostenveld
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

border = axis;

x = [ border(1) pos(1)
      border(2) pos(1) ];

y = [ pos(2)    border(3)
      pos(2)    border(4) ];

if (~ishold)
 hold on
 h = line(x, y, varargin{:});
 hold off
else
 h = line(x, y, varargin{:});
end

% make both lines the same color
set(h(2), 'Color', get(h(1), 'Color'));

