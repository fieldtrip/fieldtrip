function [varargout] = ft_plot_text(X, Y, str, varargin)

% FT_PLOT_TEXT helper function for plotting text, which can also be used in
% combination with the multiple channel layout display in FieldTrip.
%
% Use as
%   ft_plot_text(X, Y, ...)
% where optional input arguments should come in key-value pairs and may
% include
%   hpos
%   vpos
%   width
%   height
%   hlim
%   vlim
%   Color
%   FontSize
%   FontName
%   HorizontalAlignment

% Copyrights (C) 2009, Robert Oostenveld
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
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'hlim', 'vlim', 'Color', ....
  'FontSize', 'FontName', 'HorizontalAlignment','rotation','VerticalAlignment'});
hpos        = keyval('hpos',      varargin);
vpos        = keyval('vpos',      varargin);
width       = keyval('width',     varargin);
height      = keyval('height',    varargin);
hlim        = keyval('hlim',      varargin);
vlim        = keyval('vlim',      varargin);
Color       = keyval('Color',     varargin);  if isempty(Color), Color = 'k'; end
FontSize    = keyval('FontSize',  varargin);
FontName    = keyval('FontName',  varargin);
HorizontalAlignment = keyval('HorizontalAlignment',  varargin); if isempty(HorizontalAlignment), HorizontalAlignment = 'center'; end
rotation    = keyval('rotation',  varargin);  if isempty(rotation), rotation = 0; end
VerticalAlignment = keyval('VerticalAlignment',  varargin); if isempty(VerticalAlignment), VerticalAlignment = 'middle'; end

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

h = text(X, Y, str);
set(h, 'HorizontalAlignment', HorizontalAlignment);
set(h, 'Color', Color);
set(h, 'rotation', rotation);
set(h, 'VerticalAlignment',VerticalAlignment);
if ~isempty(FontSize), set(h, 'FontSize', FontSize); end
if ~isempty(FontName), set(h, 'FontName', FontName); end

% the (optional) output is the handle
if nargout == 1;
  varargout{1} = h;
end

warning(ws); %revert to original state
