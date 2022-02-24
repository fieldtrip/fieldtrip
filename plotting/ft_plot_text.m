function [varargout] = ft_plot_text(X, Y, str, varargin)

% FT_PLOT_TEXT helper function for plotting text, which can also be used in
% combination with the multiple channel layout display in FieldTrip.
%
% Use as
%   ft_plot_text(X, Y, str, ...)
%
% Optional arguments should come in key-value pairs and can include
%   'fontcolor'           = string, color specification (default = 'k')
%   'fontsize'            = number, sets the size of the text (default = 10)
%   'fontunits'           =
%   'fontname'            =
%   'fontweight'          =
%   'horizontalalignment' =
%   'verticalalignment'   =
%   'interpreter'         = string, can be 'none', 'tex' or 'latex' (default = 'none')
%   'rotation'            =
%   'tag'                 = string, the name assigned to the object. All tags with the same name can be deleted in a figure, without deleting other parts of the figure.
%
% It is possible to plot the object in a local pseudo-axis (c.f. subplot), which is specfied as follows
%   'hpos'                = horizontal position of the center of the local axes
%   'vpos'                = vertical position of the center of the local axes
%   'width'               = width of the local axes
%   'height'              = height of the local axes
%   'hlim'                = horizontal scaling limits within the local axes
%   'vlim'                = vertical scaling limits within the local axes
%
% Example
%   figure
%   ft_plot_vector(randn(1,10), rand(1,10), 'hpos', 1, 'vpos', 1, 'width', 0.2, 'height', 0.2, 'box', true)
%   ft_plot_text(0, 0 , '+',                'hpos', 1, 'vpos', 1, 'width', 0.2, 'height', 0.2)
%   axis([0 2 0 2])
%
% See also FT_PLOT_VECTOR, FT_PLOT_MATRIX, FT_PLOT_LINE, FT_PLOT_BOX

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
hpos                = ft_getopt(varargin, 'hpos');
vpos                = ft_getopt(varargin, 'vpos');
width               = ft_getopt(varargin, 'width');
height              = ft_getopt(varargin, 'height');
hlim                = ft_getopt(varargin, 'hlim');
vlim                = ft_getopt(varargin, 'vlim');
tag                 = ft_getopt(varargin, 'tag', '');
% these have to do with the font
color               = ft_getopt(varargin, 'fontcolor', 'k');
fontsize            = ft_getopt(varargin, 'fontsize',   get(0, 'defaulttextfontsize'));
fontname            = ft_getopt(varargin, 'fontname',   get(0, 'defaulttextfontname'));
fontweight          = ft_getopt(varargin, 'fontweight', get(0, 'defaulttextfontweight'));
fontunits           = ft_getopt(varargin, 'fontunits',  get(0, 'defaulttextfontunits'));
% these also have to do with the font
horizontalalignment = ft_getopt(varargin, 'horizontalalignment', 'center');
verticalalignment   = ft_getopt(varargin, 'verticalalignment', 'middle');
rotation            = ft_getopt(varargin, 'rotation', 0);
interpreter         = ft_getopt(varargin, 'interpreter', 'none');

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
    ft_warning('use hlim/vlim when specifying local axes');
    hlim = get(gca, 'XLim');
  end
  
  if isempty(vlim)
    ft_warning('use hlim/vlim when specifying local axes');
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

  X = X - hlim(1);
  X = X ./ (hlim(2)-hlim(1));
  X = X .* width;
  X = X + hpos - width/2;
  
  Y = Y - vlim(1);
  Y = Y ./ (vlim(2)-vlim(1));
  Y = Y .* height;
  Y = Y + vpos - height/2;
  
end % shortcut

% it fails on single inputs
X = double(X);
Y = double(Y);

h = text(X, Y, str, 'color', color, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
if ~isempty(horizontalalignment), set(h, 'horizontalalignment', horizontalalignment); end
if ~isempty(verticalalignment), set(h, 'verticalalignment', verticalalignment); end
if ~isempty(rotation), set(h, 'rotation', rotation); end

set(h, 'tag', tag);
set(h, 'interpreter', interpreter);

% the (optional) output is the handle
if nargout == 1
  varargout{1} = h;
end
