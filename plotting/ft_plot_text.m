function [varargout] = ft_plot_text(X, Y, str, varargin)

% FT_PLOT_TEXT helper function for plotting text, which can also be used in
% combination with the multiple channel layout display in FieldTrip.
%
% Use as
%   ft_plot_text(X, Y, str, ...)
%
% Optional arguments should come in key-value pairs and can include
%   Color               =
%   FontSize            =
%   FontName            =
%   HorizontalAlignment =
%   tag                 =
%
% It is possible to plot the object in a local pseudo-axis (c.f. subplot), which is specfied as follows
%   hpos        = horizontal position of the center of the local axes
%   vpos        = vertical position of the center of the local axes
%   width       = width of the local axes
%   height      = height of the local axes
%   hlim        = horizontal scaling limits within the local axes
%   vlim        = vertical scaling limits within the local axes

% Copyrights (C) 2009-2011, Robert Oostenveld
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

ws = warning('on', 'MATLAB:divideByZero');

% get the optional input arguments
hpos                = ft_getopt(varargin, 'hpos');
vpos                = ft_getopt(varargin, 'vpos');
width               = ft_getopt(varargin, 'width');
height              = ft_getopt(varargin, 'height');
hlim                = ft_getopt(varargin, 'hlim');
vlim                = ft_getopt(varargin, 'vlim');
Color               = ft_getopt(varargin, 'Color', 'k');
FontSize            = ft_getopt(varargin, 'FontSize');
FontName            = ft_getopt(varargin, 'FontName');
FontUnits           = ft_getopt(varargin, 'FontUnits');
HorizontalAlignment = ft_getopt(varargin, 'HorizontalAlignment', 'center');
rotation            = ft_getopt(varargin, 'rotation', 0);
VerticalAlignment   = ft_getopt(varargin, 'VerticalAlignment', 'middle');
tag                 = ft_getopt(varargin, 'tag', '');
interpreter         = ft_getopt(varargin, 'interpreter', 'tex');

if isempty(hlim) && isempty(vlim) && isempty(hpos) && isempty(vpos) && isempty(height) && isempty(width)
  % no scaling is needed, the input X and Y are already fine
  % use a shortcut to speed up the plotting
  
else
  % use the full implementation
  abc = axis;
  if isempty(hlim)
    hlim = abc([1 2]);
  end
  
  if isempty(vlim)
    vlim = abc([3 4]);
  end
  
  if isempty(hpos);
    hpos = (hlim(1)+hlim(2))/2;
  end
  
  if isempty(vpos);
    vpos = (vlim(1)+vlim(2))/2;
  end
  
  if isempty(width),
    width = hlim(2)-hlim(1);
  end
  
  if isempty(height),
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

% it fails on single inputs
X = double(X);
Y = double(Y);

h = text(X, Y, str);
set(h, 'HorizontalAlignment', HorizontalAlignment);
set(h, 'Color', Color);
set(h, 'rotation', rotation);
set(h, 'VerticalAlignment',VerticalAlignment); 
if ~isempty(FontSize),  set(h, 'FontSize', FontSize);  end
if ~isempty(FontName),  set(h, 'FontName', FontName);  end
if ~isempty(FontUnits), set(h, 'FontUnits', FontUnits); end
set(h, 'tag', tag);
set(h, 'interpreter', interpreter);

% the (optional) output is the handle
if nargout == 1;
  varargout{1} = h;
end

warning(ws); %revert to original state
