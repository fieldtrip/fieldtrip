function [h] = crosshair(pos, varargin)

% This function is DEPRECATED, please use FT_PLOT_CROSSHAIR instead.

% Copyright (C) 2003, Robert Oostenveld
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

parent = ft_getopt(varargin, 'parent');
h      = ft_getopt(varargin, 'handle');
if ~isempty(h)
  % make the parent to the first handle current axes
  set(gcf,'currentaxes',get(h(1),'parent'));
elseif ~isempty(parent)
  % make parent current axes
  set(gcf,'currentaxes',parent);
end
border = [get(gca, 'xlim') get(gca, 'ylim') get(gca, 'zlim')];

if numel(pos)==2
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
    set(h(1),'xdata',x(:,1)');
    set(h(1),'ydata',y(:,1)');
    set(h(2),'xdata',x(:,2)');
    set(h(2),'ydata',y(:,2)');
  end

  % make both lines the same color
  set(h(2), 'Color', get(h(1), 'Color'));
elseif numel(pos)==3
    
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
    set(h(1),'xdata',x(:,1)');
    set(h(1),'ydata',y(:,1)');
    set(h(1),'zdata',z(:,1)');
    set(h(2),'xdata',x(:,2)');
    set(h(2),'ydata',y(:,2)');
    set(h(2),'zdata',z(:,2)');
    set(h(3),'xdata',x(:,3)');
    set(h(3),'ydata',y(:,3)');
    set(h(3),'zdata',z(:,3)');
  end
  
  % make all lines the same color
  set(h(2), 'Color', get(h(1), 'Color'));
  set(h(3), 'Color', get(h(1), 'Color'));
end
