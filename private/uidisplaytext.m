function uidisplaytext(str, title)

% UIDISPLAYTEXT opens a figure for displaying multi-line text
% in an "edit" user interface control element.
%
% Use as
%   uidisplaytext(str, title)

% Copyright (C) 2009, Robert Oostenveld
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

if nargin<2
  title = 'unknown';
end

f = figure;
h = uicontrol('style', 'edit');

set(f, 'toolBar', 'none')
set(f, 'menuBar', 'none')

set(f, 'Name', title)
set(f, 'NumberTitle', 'off')

set(h, 'Units', 'normalized');
set(h, 'Position', [0 0 1 1]);
set(h, 'backgroundColor', [1 1 1]);
set(h, 'HorizontalAlign', 'left');
set(h, 'max', 2);
set(h, 'min', 0);
set(h, 'FontName', 'Courier');
set(h, 'FontSize', 12);
set(h, 'string', str);
