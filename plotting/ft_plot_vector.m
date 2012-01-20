function [varargout] = ft_plot_vector(varargin)

% FT_PLOT_VECTOR
%
% Use as
%   ft_plot_vector(Y, ...)
% or as
%   ft_plot_vector(X, Y, ...)
% where X and Y are similar as the input to the Matlab plot function.
%
% Optional arguments should come in key-value pairs and can include
%   axis            = draw the local axis,  can be 'yes', 'no', 'xy', 'x' or 'y'
%   box             = draw a box around the local axes, can be 'yes' or 'no'
%   highlight       = a logical vector of size Y, where 1 means that the
%                     corresponding values in Y are highlighted
%                     (according to the highlightstyle)
%   highlightstyle  = can be 'box', 'thickness', 'saturation' 
%                     ('opacity' is not supported yet, default='box')
%   tag             = a name this vector gets. All tags with the same name
%                     can be deleted in a figure, without deleting other 
%                     parts of the figure
%   color           = see MATLAB standard Line Properties
%   linewidth       = see MATLAB standard Line Properties
%   markersize      = see MATLAB standard Line Properties
%   markerfacecolor = see MATLAB standard Line Properties
%   style           = see MATLAB standard Line Properties
%   label           = see MATLAB standard Line Properties
%   fontsize        = see MATLAB standard Line Properties
%
% It is possible to plot the object in a local pseudo-axis (c.f. subplot), which is specfied as follows
%   hpos        = horizontal position of the center of the local axes
%   vpos        = vertical position of the center of the local axes
%   width       = width of the local axes
%   height      = height of the local axes
%   hlim        = horizontal scaling limits within the local axes
%   vlim        = vertical scaling limits within the local axes
%
% Example use
%   ft_plot_vector(randn(1,100), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0)

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

if nargin>1 && all(cellfun(@isnumeric, varargin(1:2)) | cellfun(@islogical, varargin(1:2)))
  % the function was called like plot(x, y, ...)
  hdat = varargin{1};
  vdat = varargin{2};
  varargin = varargin(3:end);
else
  % the function was called like plot(y, ...)
  vdat = varargin{1};
  hdat = 1:size(vdat,1);
  varargin = varargin(2:end);
end

if any(size(vdat)==1)
  % ensure that it is a column vector
  vdat = vdat(:);
end

% get the optional input arguments
hpos            = ft_getopt(varargin, 'hpos');
vpos            = ft_getopt(varargin, 'vpos');
width           = ft_getopt(varargin, 'width');
height          = ft_getopt(varargin, 'height');
hlim            = ft_getopt(varargin, 'hlim', 'maxmin');
vlim            = ft_getopt(varargin, 'vlim', 'maxmin');
style           = ft_getopt(varargin, 'style', '-');
label           = ft_getopt(varargin, 'label');
fontsize        = ft_getopt(varargin, 'fontsize');
axis            = ft_getopt(varargin, 'axis', false);
box             = ft_getopt(varargin, 'box', false);
color           = ft_getopt(varargin, 'color');
linewidth       = ft_getopt(varargin, 'linewidth', 0.5);
highlight       = ft_getopt(varargin, 'highlight');
highlightstyle  = ft_getopt(varargin, 'highlightstyle', 'box');
markersize      = ft_getopt(varargin, 'markersize', 6);
markerfacecolor = ft_getopt(varargin, 'markerfacecolor', 'none');
tag            = ft_getopt(varargin, 'tag', '');

if ~isempty(highlight) && any(size(highlight)==1)
    % ensure that it is a column vector
    highlight = highlight(:);
end
  
if ~isempty(highlight) && ~isequal(size(highlight), size(vdat))
  error('the dimensions of the highlight should be identical to the dimensions of the data');
end

% convert the yes/no strings into boolean values
box  = istrue(box);

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
      error('unsupported option for hlim')
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
      error('unsupported option for vlim')
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


% plotting lines
if isempty(color)
  h = plot(hdat, vdat, style, 'LineWidth', linewidth,'markersize',markersize,'markerfacecolor',markerfacecolor);
else
  h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'Color', color,'markersize',markersize,'markerfacecolor',markerfacecolor);
end


if ~isempty(highlight)
  if ~islogical(highlight)
    if ~all(highlight==0 | highlight==1)
      % only warn if really different from 0/1
      warning('converting mask to logical values')
    end
    highlight=logical(highlight);
  end
  
  switch highlightstyle
    case 'box'
      % find the sample number where the highlight begins and ends
      begsample = find(diff([0;highlight;0])== 1);
      endsample = find(diff([0;highlight;0])==-1)-1;
      for i=1:length(begsample)
        begx = hdat(begsample(i));
        endx = hdat(endsample(i));
        ft_plot_box([begx endx vpos-height/2 vpos+height/2], 'facecolor', [.6 .6 .6], 'edgecolor', 'none');
      end
      % plotting lines again, otherwise box will always be on top
      if isempty(color)
        h = plot(hdat, vdat, style, 'LineWidth', linewidth,'markersize',markersize,'markerfacecolor',markerfacecolor);
      else
        h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'Color', color,'markersize',markersize,'markerfacecolor',markerfacecolor);
      end
    case 'thickness'
      % find the sample number where the highligh begins and ends
      begsample = find(diff([0;highlight;0])== 1);
      endsample = find(diff([0;highlight;0])==-1)-1;
      linecolor = get(h,'Color'); % get current line color
      for i=1:length(begsample)
        hor = hdat(begsample(i):endsample(i));
        ver = vdat(begsample(i):endsample(i));
        plot(hor,ver,'linewidth',4*linewidth,'linestyle','-','Color', linecolor); % changed 3* to 4*, as 3* appeared to have no effect
      end
    case 'saturation'
      % find the sample number where the highligh begins and ends
      highlight = ~highlight; % invert the mask
      begsample = find(diff([0;highlight;0])== 1);
      endsample = find(diff([0;highlight;0])==-1)-1;
      linecolor = get(h,'Color'); % get current line color
      linecolor = (linecolor * 0.2) + 0.8; % change saturation of color
      for i=1:length(begsample)
        hor = hdat(begsample(i):endsample(i));
        ver = vdat(begsample(i):endsample(i));
        plot(hor,ver,'color',linecolor);
      end
    case 'opacity'
      error('unsupported highlightstyle')
    otherwise
      error('unsupported highlightstyle')
  end % switch highlightstyle
end


if ~isempty(label)
  boxposition(1) = hpos - width/2;
  boxposition(2) = hpos + width/2;
  boxposition(3) = vpos - height/2;
  boxposition(4) = vpos + height/2;
  h = text(boxposition(1), boxposition(4), label);
  if ~isempty(fontsize)
    set(h, 'Fontsize', fontsize);
  end
end

if box
  boxposition = zeros(1,4);
  % this plots a box around the original hpos/vpos with appropriate width/height
  x1 = hpos - width/2;
  x2 = hpos + width/2;
  y1 = vpos - height/2;
  y2 = vpos + height/2;
  
  X = [x1 x2 x2 x1 x1];
  Y = [y1 y1 y2 y2 y1];
  line(X, Y);
  
  %   % this plots a box around the original hpos/vpos with appropriate width/height
  %   boxposition(1) = hpos - width/2;
  %   boxposition(2) = hpos + width/2;
  %   boxposition(3) = vpos - height/2;
  %   boxposition(4) = vpos + height/2;
  %   ft_plot_box(boxposition, 'facecolor', 'none', 'edgecolor', 'k');
  
  % this plots a box around the complete data
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
      error('invalid specification of the "axis" option')
  end
  
  % determine where the original [0, 0] in the data is located in the scaled and shifted axes
  x0 = interp1(hlim, hpos + [-width/2  width/2 ], 0, 'linear', 'extrap');
  y0 = interp1(vlim, vpos + [-height/2 height/2], 0, 'linear', 'extrap');
  
  if xaxis
    X = [hpos-width/2  hpos+width/2];
    Y = [y0 y0];
    ft_plot_line(X, Y);
    % str = sprintf('%g', hlim(1)); ft_plot_text(X(1), Y(1), str);
    % str = sprintf('%g', hlim(2)); ft_plot_text(X(2), Y(2), str);
  end
  
  if yaxis
    X = [x0 x0];
    Y = [vpos-height/2 vpos+height/2];
    ft_plot_line(X, Y);
    % str = sprintf('%g', vlim(1)); ft_plot_text(X(1), Y(1), str);
    % str = sprintf('%g', vlim(2)); ft_plot_text(X(2), Y(2), str);
  end
end

set(h,'tag',tag);

% the (optional) output is the handle
if nargout == 1;
  varargout{1} = h;
end

if ~holdflag
  hold off
end

warning(ws); %revert to original state
