function ft_plot_crosshair(pos, varargin)

% FT_PLOT_CROSSHAIR plots a crosshair in two or three dimensions
%
% Use as
%   ft_plot_crosshair(pos, ...)
% where pos is the desired position of the crosshair.
%
% Optional input arguments should be specified in key-value pairs and can include
%   'color'    = [r g b] value or string, see PLOT
%
% Example
%   ft_plot_crosshair([0.5 0.5], 'color', 'r')

% Copyright (C) 2014, Robert Oostenveld
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

ax  = axis;

minx = ax(1);
maxx = ax(2);
miny = ax(3);
maxy = ax(4);

switch numel(pos)
  case 2
    X = [minx maxx];
    Y = pos([2 2]);
    line(X, Y, varargin{:});
    
    X = pos([1 1]);
    Y = [miny maxy];
    line(X, Y, varargin{:});
    
  case 3
    if numel(axis)==4
      % the figure is in 2D mode
      [az, el] = view;
      view([1 1 1]);
      ax = axis;
      view(az, el);
    end
    
    minz = ax(5);
    maxz = ax(6);
    
    X = [minx maxx];
    Y = pos([2 2]);
    Z = pos([3 3]);
    line(X, Y, Z, varargin{:});
    
    X = pos([1 1]);
    Y = [miny maxy];
    Z = pos([3 3]);
    line(X, Y, Z, varargin{:});
    
    X = pos([1 1]);
    Y = pos([2 2]);
    Z = [minz maxz];
    line(X, Y, Z, varargin{:});
    
end
