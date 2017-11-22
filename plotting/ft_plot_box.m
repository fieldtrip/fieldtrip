function [varargout] = ft_plot_box(position, varargin)

% FT_PLOT_BOX plots the outline of a box that is specified by its lower
% left and upper right corner
%
% Use as
%   ft_plot_box(position, ...)
% where the position of the box is specified as is [x1, x2, y1, y2].
%
% Optional arguments should come in key-value pairs and can include
%   'facealpha'       = transparency value between 0 and 1
%   'facecolor'       = color specification as [r g b] values or a string, for example 'brain', 'cortex', 'skin', 'red', 'r'
%   'edgecolor'       = color specification as [r g b] values or a string, for example 'brain', 'cortex', 'skin', 'red', 'r'
%   'tag'             = string, the name assigned to the object. All tags with the same name can be deleted in a figure, without deleting other parts of the figure.
%
% It is possible to plot the object in a local pseudo-axis (c.f. subplot), which is specfied as follows
%   'hpos'            = horizontal position of the center of the local axes
%   'vpos'            = vertical position of the center of the local axes
%   'width'           = width of the local axes
%   'height'          = height of the local axes
%   'hlim'            = horizontal scaling limits within the local axes
%   'vlim'            = vertical scaling limits within the local axes
%   'parent'          = handle which is set as the parent for all plots
%
% Example
%   ft_plot_box([-1 1 2 3], 'facecolor', 'b')
%   axis([-4 4 -4 4])
%
% See also FT_PLOT_LINE

% Copyrights (C) 2009-2011, Robert Oostenveld
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

ws = warning('on', 'MATLAB:divideByZero');

% get the optional input arguments
hpos        = ft_getopt(varargin, 'hpos');
vpos        = ft_getopt(varargin, 'vpos');
width       = ft_getopt(varargin, 'width');
height      = ft_getopt(varargin, 'height');
hlim        = ft_getopt(varargin, 'hlim');
vlim        = ft_getopt(varargin, 'vlim');
facealpha   = ft_getopt(varargin, 'facealpha', 1);
facecolor   = ft_getopt(varargin, 'facecolor', 'none');
edgecolor   = ft_getopt(varargin, 'edgecolor', 'k');
tag         = ft_getopt(varargin, 'tag',       '');
parent      = ft_getopt(varargin, 'parent', []);

% color management
if ischar(facecolor) && exist([facecolor '.m'], 'file')
	facecolor = eval(facecolor);
end
if ischar(edgecolor) && exist([edgecolor '.m'], 'file')
	edgecolor = eval(edgecolor);
end

% convert the two cornerpoints into something that the patch function understands
% the box position is represented just like the argument to the AXIS function
x1 = position(1);
x2 = position(2);
y1 = position(3);
y2 = position(4);
X = [x1 x2 x2 x1 x1];
Y = [y1 y1 y2 y2 y1];

if isempty(hlim) && isempty(vlim) && isempty(hpos) && isempty(vpos) && isempty(height) && isempty(width)
  % no scaling is needed, the input X and Y are already fine
  % use a shortcut to speed up the plotting
  
else
  % use the full implementation
  if isempty(hlim)
    hlim = get(gca, 'XLim');
  end
  
  if isempty(vlim)
    vlim = get(gca, 'YLim');
  end
  
  if isempty(hpos)
    hpos = (hlim(1)+hlim(2))/2;
  end
  
  if isempty(vpos)
    vpos = (vlim(1)+vlim(2))/2;
  end
  
  if isempty(width)
    width = hlim(2)-hlim(1);
  end
  
  if isempty(height)
    height = vlim(2)-vlim(1);
  end
  
  % first shift the horizontal axis to zero
  X = X - (hlim(1)+hlim(2))/2;
  % then scale to length 1
  X = X ./ (hlim(2)-hlim(1));
  % then scale to the new width
  X = X .* width;
  % then shift to the new horizontal position
  X = X + hpos;
  
  % first shift the vertical axis to zero
  Y = Y - (vlim(1)+vlim(2))/2;
  % then scale to length 1
  Y = Y ./ (vlim(2)-vlim(1));
  % then scale to the new width
  Y = Y .* height;
  % then shift to the new vertical position
  Y = Y + vpos;
  
end % shortcut

% use an arbitrary color, which will be replaced by the correct color a few lines down
C = 0;

h = patch(X, Y, C);
set(h, 'FaceAlpha', facealpha)
set(h, 'FaceColor', facecolor)
set(h, 'EdgeColor', edgecolor)
set(h, 'tag', tag);

if ~isempty(parent)
  set(h, 'Parent', parent);
end

% the (optional) output is the handle
if nargout == 1
  varargout{1} = h;
end

warning(ws); % revert to original state

