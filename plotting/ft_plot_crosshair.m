function h = ft_plot_crosshair(pos, varargin)

% FT_PLOT_CROSSHAIR plots a crosshair at a specified position in two [x, y] or three
% [x, y, z] dimensions. 
%
% Use as
%   h = ft_plot_crosshair(pos, ...)
% where pos is the desired position of the crosshair. The handles of the lines are
% returned.
%
% Optional input arguments should be specified in key-value pairs and can include
%   'color'    = [r g b] value or string, see PLOT
%   'parent'   = handle of the parent axes
%   'handle'   = handle of the existing line objects to be updated
% 
% You can specify the handles of existing line objects which will be then updated,
% rather than creating a new set of lines. If both parent and handle ar specified,
% the handle option prevail.
% 
% Example
%   ft_plot_crosshair([0.5 0.5], 'color', 'r')
%
% See also TEXT, LINE

% Copyright (C) 2003-2017, Robert Oostenveld
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

% get the optional input arguments
color  = ft_getopt(varargin, 'color');
parent = ft_getopt(varargin, 'parent');
h      = ft_getopt(varargin, 'handle');

if ~isempty(h)
  % the parent to the first handle is set as the current axes
  set(gcf, 'currentaxes', get(h(1),'parent'));
elseif ~isempty(parent)
  % the parent is set as the current axes
  set(gcf, 'currentaxes', parent);
else
  % the current axes stay as they are
end

% color management
if ~isempty(color)
  if ischar(color) && exist([color '.m'], 'file')
    color = eval(color);
  end
  varargin = ft_setopt(varargin, 'color', color);
end

% determine the size of the figure
border = [get(gca, 'xlim') get(gca, 'ylim') get(gca, 'zlim')];

switch numel(pos)
  case 2
    x = [ border(1) pos(1)
          border(2) pos(1) ];
    y = [ pos(2)    border(3)
          pos(2)    border(4) ];
    
    if isempty(h) && (~ishold)
      hold on
      h = line(x, y, varargin{:});
      hold off
    elseif isempty(h)
      h = line(x, y, varargin{:});
    else
      set(h(1), 'xdata', x(:,1)');
      set(h(1), 'ydata', y(:,1)');
      set(h(2), 'xdata', x(:,2)');
      set(h(2), 'ydata', y(:,2)');
    end
    
    if ~isempty(color)
      set(h(1), 'color', color);
      set(h(2), 'color', color);
    else
      % ensure they have the same color
      set(h(2), 'color', get(h(1), 'color'));
    end
    
  case 3
    x = [ border(1) pos(1)    pos(1)
          border(2) pos(1)    pos(1)];
    y = [ pos(2)    border(3) pos(2)
          pos(2)    border(4) pos(2)];
    z = [ pos(3)    pos(3)    border(5)
          pos(3)    pos(3)    border(6)];
    
    if isempty(h) && (~ishold)
      hold on
      h = line(x, y, z, varargin{:});
      hold off
    elseif isempty(h)
      h = line(x, y, z, varargin{:});
    else
      set(h(1), 'xdata', x(:,1)');
      set(h(1), 'ydata', y(:,1)');
      set(h(1), 'zdata', z(:,1)');
      set(h(2), 'xdata', x(:,2)');
      set(h(2), 'ydata', y(:,2)');
      set(h(2), 'zdata', z(:,2)');
      set(h(3), 'xdata', x(:,3)');
      set(h(3), 'ydata', y(:,3)');
      set(h(3), 'zdata', z(:,3)');
    end
    
    if ~isempty(color)
      set(h(1), 'color', color);
      set(h(2), 'color', color);
      set(h(3), 'color', color);
    else
      % ensure they have the same color
      set(h(2), 'color', get(h(1), 'color'));
      set(h(3), 'color', get(h(1), 'color'));
    end
    
end % switch
