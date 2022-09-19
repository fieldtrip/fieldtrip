function h = ft_plot_line(X, Y, varargin)

% FT_PLOT_LINE helper function for plotting a line, which can also be used in
% combination with the multiple channel layout display in FieldTrip.
%
% Use as
%   ft_plot_line(X, Y, ...)
%
% Optional arguments should come in key-value pairs and can include
%   'color'           =
%   'linestyle'       =
%   'linewidth'       =
%   'tag'             = string, the name assigned to the object. All tags with the same name can be deleted in a figure, without deleting other parts of the figure.
%
% It is possible to plot the object in a local pseudo-axis (c.f. subplot), which is specfied as follows
%   'hpos'            = horizontal position of the center of the local axes
%   'vpos'            = vertical position of the center of the local axes
%   'width'           = width of the local axes
%   'height'          = height of the local axes
%   'hlim'            = horizontal scaling limits within the local axes
%   'vlim'            = vertical scaling limits within the local axes
%
% See also FT_PLOT_BOX, FT_PLOT_CROSSHAIR

% Copyrights (C) 2009-2022, Robert Oostenveld
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
hpos        = ft_getopt(varargin, 'hpos');
vpos        = ft_getopt(varargin, 'vpos');
width       = ft_getopt(varargin, 'width');
height      = ft_getopt(varargin, 'height');
hlim        = ft_getopt(varargin, 'hlim');
vlim        = ft_getopt(varargin, 'vlim');
color       = ft_getopt(varargin, 'color',      'k');
linestyle   = ft_getopt(varargin, 'linestyle',  '-');
linewidth   = ft_getopt(varargin, 'linewidth',  0.5);
tag         = ft_getopt(varargin, 'tag',        '');

% color management
if ischar(color) && exist([color '.m'], 'file')
  color = eval(color);
elseif ischar(color) && ismember(color, htmlcolors)
  color = htmlcolors(color);
end

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

h = line(X, Y, 'Color', color, 'LineStyle', linestyle, 'LineWidth', linewidth);
set(h, 'tag', tag);
