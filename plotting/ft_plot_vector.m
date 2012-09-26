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
%   highlight       = a logical vector of size Y, where 1 means that the corresponding values in Y are highlighted (according to the highlightstyle)
%   highlightstyle  = can be 'box', 'thickness', 'saturation', 'difference' (default='box')
%   tag             = a name this vector gets. All tags with the same name can be deleted in a figure, without deleting other parts of the figure
%   color           = see MATLAB standard Line Properties, it can also be for example 'rbg' to plot three lines with different colors
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

% FIXME if the color is Nx3, it plots the line segments with gradual
% colors, this needs to be cleaned up and improved. At the moment it
% happens around line 324.
%
% Example
%   colormap hot; rgb = colormap; rgb = interp1(1:64, rgb, linspace(1,64,100));
%   ft_plot_vector(1:100, 'color', rgb);

% Copyrights (C) 2009-2012, Robert Oostenveld
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
  if iscolumn(vdat),
    hdat = 1:size(vdat,1);
  else
    hdat = 1:size(vdat,2);
  end
  varargin = varargin(2:end);
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
tag             = ft_getopt(varargin, 'tag', '');
parent          = ft_getopt(varargin, 'parent', []);

% if any(size(vdat)==1)
%   % ensure that it is a column vector
%   vdat = vdat(:);
% end

% if ~isempty(highlight) && any(size(highlight)==1)
%     % ensure that it is a column vector
%     highlight = highlight(:);
% end

if ischar(color)
  % it should be a column array
  color = color(:);
end

if ~isempty(highlight) && size(highlight,2)~=size(vdat,2)
  error('the dimensions of the highlight should be identical to the dimensions of the data');
end

% convert the yes/no strings into boolean values
box = istrue(box);

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

if ~isempty(highlight) && ~islogical(highlight)
  if ~all(highlight==0 | highlight==1)
    % only warn if really different from 0/1
    warning('converting mask to logical values')
  end
  highlight=logical(highlight);
end

switch highlightstyle
  case 'box'
    % find the sample number where the highlight begins and ends
    begsample = find(diff([0 highlight 0])== 1);
    endsample = find(diff([0 highlight 0])==-1)-1;
    for i=1:length(begsample)
      begx = hdat(begsample(i));
      endx = hdat(endsample(i));
      ft_plot_box([begx endx vpos-height/2 vpos+height/2], 'facecolor', [.6 .6 .6], 'edgecolor', 'none', 'parent', parent);
    end
    
  case 'thickness'
    % find the sample number where the highligh begins and ends
    begsample = find(diff([0 highlight 0])== 1);
    endsample = find(diff([0 highlight 0])==-1)-1;
    for j=1:size(vdat,1)
      for i=1:length(begsample)
        hor = hdat(   begsample(i):endsample(i));
        ver = vdat(j, begsample(i):endsample(i));
        if isempty(color)
          plot(hor,ver,'linewidth',4*linewidth,'linestyle','-'); % changed 3* to 4*, as 3* appeared to have no effect
        elseif ischar(color) && size(color,1)==1
          % plot all lines with the same color
          plot(hor,ver,'linewidth',4*linewidth,'linestyle','-','Color', color); % changed 3* to 4*, as 3* appeared to have no effect
        else
          % plot each line with its own color
          plot(hor,ver,'linewidth',4*linewidth,'linestyle','-','Color', color(j)); % changed 3* to 4*, as 3* appeared to have no effect
        end
        
      end
    end
    
  case 'saturation'
    % find the sample number where the highligh begins and ends
    highlight = ~highlight; % invert the mask
    begsample = find(diff([0 highlight 0])== 1);
    endsample = find(diff([0 highlight 0])==-1)-1;
    % start with plotting the lines
    for i=1:size(vdat,1)
      if isempty(color)
        h = plot(hdat, vdat, style, 'LineWidth', linewidth,'markersize',markersize,'markerfacecolor',markerfacecolor);
      elseif ischar(color) && size(color,1)==1
        % plot all lines with the same color
        h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'Color', color, 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      else
        % plot each line with its own color
        h = plot(hdat, vdat(i,:), style, 'LineWidth', linewidth, 'Color', color(i), 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      end
      if ~isempty(parent)
        set(h, 'Parent', parent);
      end
      linecolor = get(h, 'color');
      linecolor = (linecolor * 0.2) + 0.8; % change saturation of color
      for j=1:length(begsample)
        hor = hdat(   begsample(j):endsample(j));
        ver = vdat(i, begsample(j):endsample(j));
        h = plot(hor,ver,'color',linecolor);
        
        if ~isempty(parent)
          set(h, 'Parent', parent);
        end
      end
    end
    
  case 'difference'
    if size(vdat,1)~=2
      error('this only works if exactly two lines are plotted');
    end
    hdatbeg = [hdat(:,1) (hdat(:,1:end-1) + hdat(:,2:end))/2            ];
    hdatend = [          (hdat(:,1:end-1) + hdat(:,2:end))/2 hdat(:,end)];
    vdatbeg = [vdat(:,1) (vdat(:,1:end-1) + vdat(:,2:end))/2            ];
    vdatend = [          (vdat(:,1:end-1) + vdat(:,2:end))/2 vdat(:,end)];
    begsample = find(diff([0 highlight 0])== 1);
    endsample = find(diff([0 highlight 0])==-1)-1;
    for i=1:length(begsample)
      X = [hdatbeg(1,begsample(i)) hdat(1,begsample(i):endsample(i)) hdatend(1,endsample(i)) hdatend(1,endsample(i)) hdat(1,endsample(i):-1:begsample(i)) hdatbeg(1,begsample(i))];
      Y = [vdatbeg(1,begsample(i)) vdat(1,begsample(i):endsample(i)) vdatend(1,endsample(i)) vdatend(2,endsample(i)) vdat(2,endsample(i):-1:begsample(i)) vdatbeg(2,begsample(i))];
      h = patch(X, Y, [.6 .6 .6]);
      set(h, 'linestyle', 'no');
      if ~isempty(parent)
        set(h, 'Parent', parent);
      end
    end
    
  otherwise
    % no hightlighting needs to be done
end % switch highlightstyle


switch highlightstyle
  case 'saturation'
    % this plots the lines together with the hightlights, nothing left to do
    
  otherwise
    % plot the actual lines after the highlight box or patch, otherwise those will be on top
    if isempty(color)
      h = plot(hdat, vdat, style, 'LineWidth', linewidth,'markersize',markersize,'markerfacecolor',markerfacecolor);
    elseif ischar(color) && length(color)>1
      % plot each line with its own color
      for i=1:size(vdat,1)
        h = plot(hdat, vdat(i,:), style, 'LineWidth', linewidth, 'Color', color(i), 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      end
    elseif isnumeric(color) && length(color)==length(vdat)
      for i=1:(length(vdat)-1)
        h = plot(hdat(i:i+1), vdat(i:i+1), style, 'LineWidth', linewidth, 'Color', mean(color([i i+1],:),1), 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      end
    end
     if ~isempty(parent)
       set(h, 'Parent', parent);
     end
end % switch highlightstyle


if ~isempty(label)
  boxposition(1) = hpos - width/2;
  boxposition(2) = hpos + width/2;
  boxposition(3) = vpos - height/2;
  boxposition(4) = vpos + height/2;
  h = text(boxposition(1), boxposition(4), label);
  if ~isempty(fontsize)
    set(h, 'Fontsize', fontsize);
  end
  
  if ~isempty(parent)
    set(h, 'Parent', parent);
  end
end

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
