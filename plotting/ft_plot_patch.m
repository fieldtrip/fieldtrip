function [varargout] = ft_plot_patch(hdat, vdat, varargin)

% FT_PLOT_PATCH plot a colored shape, similar to the MATLAB patch() function. It is 
% similar in usage as ft_plot_vector, and they can be combined, for example,
% to plot an area equivalent to a SEM or STD-DEV around a line.
%
% Use as
%   ft_plot_patch(X, Y, ...)
% where X and Y are similar as the input to the MATLAB patch() function.
%
% Optional arguments should come in key-value pairs and can include
%   'axis'            = draw the local axis,  can be 'yes', 'no', 'xy', 'x' or 'y'
%   'box'             = draw a box around the local axes, can be 'yes' or 'no'
%   'tag'             = string, the name assigned to the object. All tags with the same name can be deleted in a figure, without deleting other parts of the figure.
%   'facecolor'       = see MATLAB standard patch properties 
%   'facealpha'       = see MATLAB standard patch properties (note, approx. transparency can be achieved using 'facecolor')
%   'edgecolor'       = see MATLAB standard patch properties (default is 'none') (equivalent to 'linecolor' in PLOT)
%   'linestyle'       = see MATLAB standard patch properties 
%   'linewidth'       = see MATLAB standard patch properties 
%
% The color of the patchand the edges (i.e. border lines) can be specified in a variety of ways
%   - as a string with one character per line that you want to plot. Supported colors are the same as in PATCH, i.e. 'bgrcmykw'.
%   - as an 'RGB triplet', a 1x3 vector with values between 0 and 1
%   - as 'none' if you do not want the face of the patch to be filled (useful when you want to plot an empty box).
%
% It is possible to plot the object in a local pseudo-axis (c.f. subplot), which is specfied as follows
%   'hpos'            = horizontal position of the center of the local axes
%   'vpos'            = vertical position of the center of the local axes
%   'width'           = width of the local axes
%   'height'          = height of the local axes
%   'hlim'            = horizontal scaling limits within the local axes
%   'vlim'            = vertical scaling limits within the local axes
%
% Example
%   hdat = [1:10 10:-1:1];
%   vdat = rand(1,10);
%   vdat = [vdat vdat(end:-1:1)+1];
%   ft_plot_patch(hdat, vdat)
%
% See also FT_PLOT_VECTOR, PATCH, PLOT

% Copyrights (C) 2015, Roemer van der Meij
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
hpos            = ft_getopt(varargin, 'hpos');
vpos            = ft_getopt(varargin, 'vpos');
width           = ft_getopt(varargin, 'width');
height          = ft_getopt(varargin, 'height');
hlim            = ft_getopt(varargin, 'hlim', 'maxmin');
vlim            = ft_getopt(varargin, 'vlim', 'maxmin');
axis            = ft_getopt(varargin, 'axis', false);
box             = ft_getopt(varargin, 'box', false);
tag             = ft_getopt(varargin, 'tag', '');
parent          = ft_getopt(varargin, 'parent', []);
facecolor       = ft_getopt(varargin, 'facecolor', 'b');
facealpha       = ft_getopt(varargin, 'facealpha', 1);
edgecolor       = ft_getopt(varargin, 'edgecolor', 'none');
linestyle       = ft_getopt(varargin, 'linestyle', 'none');
linewidth       = ft_getopt(varargin, 'linewidth', .5);

% convert the yes/no strings into boolean values
box = istrue(box);

% color management
if ischar(facecolor) && exist([facecolor '.m'], 'file')
  facecolor = eval(facecolor);
end
if ischar(edgecolor) && exist([edgecolor '.m'], 'file')
  edgecolor = eval(edgecolor);
end

% this should be a string, because valid options include yes, no, xy, x, y
if isequal(axis, true)
  axis = 'yes';
elseif isequal(axis, false)
  axis = 'no';
end

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

if ischar(hlim)
  switch hlim
    case 'maxmin'
      hlim = [min(hdat) max(hdat)];
    case 'maxabs'
      hlim = max(abs(hdat));
      hlim = [-hlim hlim];
    otherwise
      ft_error('unsupported option for hlim')
  end % switch
end % if ischar

if ischar(vlim)
  switch vlim
    case 'maxmin'
      vlim = [min(vdat(:)) max(vdat(:))];
    case 'maxabs'
      vlim = max(abs(vdat(:)));
      vlim = [-vlim vlim];
    otherwise
      ft_error('unsupported option for vlim')
  end % switch
end % if ischar

if vlim(1)==vlim(2)
  % vertical scaling cannot be determined, behave consistent to the plot() function
  vlim = [-1 1];
end

% these must be floating point values and not integers, otherwise the scaling fails
hdat = double(hdat);
vdat = double(vdat);
hlim = double(hlim);
vlim = double(vlim);

if isempty(hpos) && ~isempty(hlim)
  hpos = (hlim(1)+hlim(2))/2;
end
if isempty(vpos) && ~isempty(vlim)
  vpos = (vlim(1)+vlim(2))/2;
end

if isempty(width) && ~isempty(hlim)
  width = hlim(2)-hlim(1);
end

if isempty(height) && ~isempty(vlim)
  height = vlim(2)-vlim(1);
end

% first shift the horizontal axis to zero
if any(hlim) ~= 0
  hdat = hdat - (hlim(1)+hlim(2))/2;
  % then scale to length 1
  if hlim(2)-hlim(1)~=0
    hdat = hdat ./ (hlim(2)-hlim(1));
  else
    hdat = hdat /hlim(1);
  end
  % then scale to the new width
  hdat = hdat .* width;
end
% then shift to the new horizontal position
hdat = hdat + hpos;

if any(vlim) ~= 0
  % first shift the vertical axis to zero
  vdat = vdat - (vlim(1)+vlim(2))/2;
  % then scale to length 1
  vdat = vdat / (vlim(2)-vlim(1));
  % then scale to the new width
  vdat = vdat .* height;
end
% then shift to the new vertical position
vdat = vdat + vpos;

% plot the patch
h = patch(hdat, vdat, facecolor, 'facealpha', facealpha, 'edgecolor', edgecolor, 'linestyle', linestyle, 'linewidth', linewidth);

if box
  % this plots a box around the original hpos/vpos with appropriate width/height
  x1 = hpos - width/2;
  x2 = hpos + width/2;
  y1 = vpos - height/2;
  y2 = vpos + height/2;
  
  X = [x1 x2 x2 x1 x1];
  Y = [y1 y1 y2 y2 y1];
  h = line(X, Y);
  set(h, 'color', 'k');
  
  % this plots a box around the original hpos/vpos with appropriate width/height
  % boxposition = zeros(1,4);
  % boxposition(1) = hpos - width/2;
  % boxposition(2) = hpos + width/2;
  % boxposition(3) = vpos - height/2;
  % boxposition(4) = vpos + height/2;
  % ft_plot_box(boxposition, 'facecolor', 'none', 'edgecolor', 'k');
  
  % this plots a box around the complete data
  % boxposition = zeros(1,4);
  % boxposition(1) = hlim(1);
  % boxposition(2) = hlim(2);
  % boxposition(3) = vlim(1);
  % boxposition(4) = vlim(2);
  % ft_plot_box(boxposition, 'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim);
end

if ~isempty(axis) && ~strcmp(axis, 'no')
  switch axis
    case {'yes' 'xy'}
      xaxis = true;
      yaxis = true;
    case {'x'}
      xaxis = true;
      yaxis = false;
    case {'y'}
      xaxis = false;
      yaxis = true;
    otherwise
      ft_error('invalid specification of the "axis" option')
  end
  
  if xaxis
    % x-axis should touch 0,0
    xrange = hlim;
    if sign(xrange(1))==sign(xrange(2))
      [dum minind] = min(abs(hlim));
      xrange(minind) = 0;
    end
    ft_plot_line(xrange, [0 0], 'hpos', hpos, 'vpos', vpos, 'hlim', hlim, 'vlim', vlim, 'width', width, 'height', height);
  end
  if yaxis
    % y-axis should touch 0,0
    yrange = vlim;
    if sign(yrange(1))==sign(yrange(2))
      [dum minind] = min(abs(vlim));
      yrange(minind) = 0;
    end
    ft_plot_line([0 0], yrange, 'hpos', hpos, 'vpos', vpos, 'hlim', hlim, 'vlim', vlim, 'width', width, 'height', height);
  end
end

set(h, 'tag', tag);

if ~isempty(parent)
  set(h, 'Parent', parent);
end

% the (optional) output is the handle
if nargout == 1;
  varargout{1} = h;
end

if ~holdflag
  hold off
end

warning(ws); % revert to original state
