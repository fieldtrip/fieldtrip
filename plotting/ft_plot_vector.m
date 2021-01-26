function [varargout] = ft_plot_vector(varargin)

% FT_PLOT_VECTOR visualizes a vector as a line, similar to PLOT.
%
% Use as
%   ft_plot_vector(Y, ...)
% or as
%   ft_plot_vector(X, Y, ...)
% where X and Y are similar as the input to the MATLAB plot function.
%
% Optional arguments should come in key-value pairs and can include
%   'axis'            = draw the local axis,  can be 'yes', 'no', 'xy', 'x' or 'y'
%   'highlight'       = a logical vector of size Y, where 1 means that the corresponding values in Y are highlighted (according to the highlightstyle)
%   'highlightstyle'  = can be 'box', 'thickness', 'saturation', 'difference' (default='box')
%   'color'           = see MATLAB standard line properties and see below
%   'facecolor'       = color for the highlighted box/difference (default = [0.6 0.6 0.6])
%   'facealpha'       = transparency for the highlighted box, between 0 and 1 (default = 1)
%   'linewidth'       = see MATLAB standard line properties
%   'markersize'      = see MATLAB standard line properties
%   'markerfacecolor' = see MATLAB standard line properties
%   'style'           = see MATLAB standard line properties
%   'tag'             = string, the name assigned to the object. All tags with the same name can be deleted in a figure, without deleting other parts of the figure.
%   'box'             = draw a box around the local axes, can be 'yes' or 'no'
%
% The line color can be specified in a variety of ways
%   - as a string with one character per line that you want to plot. Supported colors are the same as in PLOT, i.e. 'bgrcmykw'.
%   - as 'none' if you do not want the lines to be plotted (useful in combination with the difference highlightstyle).
%   - as a Nx3 matrix, where N=length(x), to use graded RGB colors along the line
%
% It is possible to plot the object in a local pseudo-axis (c.f. subplot), which is specfied as follows
%   'hpos'            = horizontal position of the center of the local axes
%   'vpos'            = vertical position of the center of the local axes
%   'width'           = width of the local axes
%   'height'          = height of the local axes
%   'hlim'            = horizontal scaling limits within the local axes
%   'vlim'            = vertical scaling limits within the local axes
%
% When using a local pseudo-axis, you can plot a label next to the data
%   'label'           = string, label to be plotted in the corner of the box
%   'labelpos'        = string, position for the label (default = 'upperleft')
%   'fontcolor'       = string, color specification (default = 'k')
%   'fontsize'        = number, sets the size of the text (default = 10)
%   'fontunits'       =
%   'fontname'        =
%   'fontweight'      =
%
% Example 1
%   subplot(2,1,1); ft_plot_vector(1:100, randn(1,100), 'color', 'r')
%   subplot(2,1,2); ft_plot_vector(1:100, randn(1,100), 'color', rand(100,3))
%
% Example 2
%   ft_plot_vector(randn(1,100), 'width', 0.9, 'height', 0.9, 'hpos', 0, 'vpos', 0, 'box', 'yes')
%   ft_plot_vector(randn(1,100), 'width', 0.9, 'height', 0.9, 'hpos', 1, 'vpos', 0, 'box', 'yes')
%   ft_plot_vector(randn(1,100), 'width', 0.9, 'height', 0.9, 'hpos', 0, 'vpos', 1, 'box', 'yes')
%
% Example 3
%  x = 1:100; y = hann(100)';
%  subplot(3,1,1); ft_plot_vector(x, y, 'highlight', y>0.8, 'highlightstyle', 'box');
%  subplot(3,1,2); ft_plot_vector(x, y, 'highlight', y>0.8, 'highlightstyle', 'thickness');
%  subplot(3,1,3); ft_plot_vector(x, y, 'highlight', y>0.8, 'highlightstyle', 'saturation');
%
% Example 4
%  x = 1:100; y = hann(100)'; ymin = 0.8*y; ymax = 1.2*y;
%  ft_plot_vector(x, [ymin; ymax], 'highlight', ones(size(y)), 'highlightstyle', 'difference', 'color', 'none');
%  ft_plot_vector(x, y);
%
% Example 5
%  colormap hot;
%  rgb = colormap;
%  rgb = interp1(1:64, rgb, linspace(1,64,100));
%  ft_plot_vector(1:100, 'color', rgb);
%
% See also FT_PLOT_MATRIX, PLOT

% Copyrights (C) 2009-2013, Robert Oostenveld
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

ws = ft_warning('on', 'MATLAB:divideByZero');

if nargin>1 && all(cellfun(@isnumeric, varargin(1:2)) | cellfun(@islogical, varargin(1:2)))
  % the function was called like plot(x, y, ...)
  hdat = varargin{1};
  vdat = varargin{2};
  varargin = varargin(3:end);
else
  % the function was called like plot(y, ...)
  vdat = varargin{1};
  if size(vdat, 1) > 1
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
labelpos        = ft_getopt(varargin, 'labelpos', 'upperleft');
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
% these have to do with the font of the label
fontcolor       = ft_getopt(varargin, 'fontcolor', 'k'); % default is black
fontsize        = ft_getopt(varargin, 'fontsize',   get(0, 'defaulttextfontsize'));
fontname        = ft_getopt(varargin, 'fontname',   get(0, 'defaulttextfontname'));
fontweight      = ft_getopt(varargin, 'fontweight', get(0, 'defaulttextfontweight'));
fontunits       = ft_getopt(varargin, 'fontunits',  get(0, 'defaulttextfontunits'));
facecolor       = ft_getopt(varargin, 'facecolor', [0.6 0.6 0.6]);
facealpha       = ft_getopt(varargin, 'facealpha', 1);

% if any(size(vdat)==1)
%   % ensure that it is a column vector
%   vdat = vdat(:);
% end

% if ~isempty(highlight) && any(size(highlight)==1)
%     % ensure that it is a column vector
%     highlight = highlight(:);
% end

npos  = numel(hdat);
nline = numel(vdat)/npos;

if ~isequal(size(hdat), [1 npos])
  hdat = hdat';
end
if ~isequal(size(vdat), [nline npos])
  vdat = vdat';
end

if ischar(color) && ~strcmp(color, 'none')
  % it should be a column array
  color = color(:);
  
  if numel(color) > nline
    % more colors specified than lines, just take the first nline
    color = color(1:nline);
  end
end

if strcmp(highlightstyle, 'difference') && isempty(highlight)
  ft_warning('highlight is empty, highlighting the whole time interval');
  highlight = ones(size(hdat));
end

if ~isempty(highlight)
  if numel(highlight)~=npos
    ft_error('the length of the highlight vector should correspond to the length of the data');
  else
    % make sure the vector points in the same direction as the data
    highlight = reshape(highlight, size(hdat));
  end
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
if any(hlim)
  hdat = hdat - (hlim(1)+hlim(2))/2;
  % then scale to length 1
  if (hlim(2)-hlim(1))~=0
    hdat = hdat ./ (hlim(2)-hlim(1));
  else
    hdat = hdat / hlim(1);
  end
  % then scale to the new width
  hdat = hdat .* width;
end
% then shift to the new horizontal position
hdat = hdat + hpos;

if any(vlim)
  % first shift the vertical axis to zero
  vdat = vdat - (vlim(1)+vlim(2))/2;
  % then scale to length 1
  vdat = vdat / (vlim(2)-vlim(1));
  % then scale to the new height
  vdat = vdat .* height;
end
% then shift to the new vertical position
vdat = vdat + vpos;

if ~isempty(highlight) && ~islogical(highlight)
  if ~all(highlight==0 | highlight==1)
    % only warn if really different from 0/1
    ft_warning('converting mask to logical values')
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
      ft_plot_box([begx endx vpos-height/2 vpos+height/2], 'facecolor', facecolor, 'facealpha', facealpha, 'edgecolor', 'none', 'parent', parent);
    end
    
  case 'thickness'
    % find the sample number where the highligh begins and ends
    begsample = find(diff([0 highlight 0])== 1);
    endsample = find(diff([0 highlight 0])==-1)-1;
    for j=1:nline
      for i=1:length(begsample)
        hor = hdat(   begsample(i):endsample(i));
        ver = vdat(j, begsample(i):endsample(i));
        if isempty(color)
          plot(hor,ver,'linewidth',4*linewidth,'linestyle','-'); % changed 3* to 4*, as 3* appeared to have no effect
        elseif ischar(color) && numel(color)==1
          % plot all lines with the same color
          plot(hor,ver,'linewidth',4*linewidth,'linestyle','-','Color', color); % changed 3* to 4*, as 3* appeared to have no effect
        elseif isnumeric(color) && isequal(size(color), [1 3])
          % plot all lines with the same RGB color
          plot(hor,ver,'linewidth',4*linewidth,'linestyle','-','Color', color); % changed 3* to 4*, as 3* appeared to have no effect
        else
          % plot each line with its own color
          plot(hor,ver,'linewidth',4*linewidth,'linestyle','-','Color', color(j)); % changed 3* to 4*, as 3* appeared to have no effect
        end
      end
    end
    
  case 'saturation'
    % find the sample number where the highlight begins and ends
    highlight = ~highlight; % invert the mask
    begsample = find(diff([0 highlight 0])== 1);
    endsample = find(diff([0 highlight 0])==-1)-1;
    % start with plotting the lines
    for i=1:nline
      if isempty(color)
        h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      elseif ischar(color) && numel(color)==1
        % plot all lines with the same color
        h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'Color', color, 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      elseif isnumeric(color) && isequal(size(color), [1 3])
        % plot all lines with the same RGB color
        h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'Color', color, 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      else
        % plot each line with its own color
        h = plot(hdat, vdat(i,:), style, 'LineWidth', linewidth, 'Color', color(i), 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      end
      linecolor = get(h, 'color');
      linecolor = (linecolor * 0.2) + 0.8; % change saturation of color
      for j=1:length(begsample)
        hor = hdat(   begsample(j):endsample(j));
        ver = vdat(i, begsample(j):endsample(j));
        h = plot(hor,ver,'color',linecolor);
      end
    end
    
  case 'difference'
    if nline~=2
      ft_error('this only works if exactly two lines are plotted');
    end
    hdatbeg = [hdat(:,1) (hdat(:,1:end-1) + hdat(:,2:end))/2            ];
    hdatend = [          (hdat(:,1:end-1) + hdat(:,2:end))/2 hdat(:,end)];
    vdatbeg = [vdat(:,1) (vdat(:,1:end-1) + vdat(:,2:end))/2            ];
    vdatend = [          (vdat(:,1:end-1) + vdat(:,2:end))/2 vdat(:,end)];
    begsample = find(diff([0 highlight 0])== 1)  ;
    endsample = find(diff([0 highlight 0])==-1)-1;
    for i=1:length(begsample)
      X = [hdatbeg(1,begsample(i)) hdat(1,begsample(i):endsample(i)) hdatend(1,endsample(i)) hdatend(1,endsample(i)) hdat(1,endsample(i):-1:begsample(i)) hdatbeg(1,begsample(i))];
      Y = [vdatbeg(1,begsample(i)) vdat(1,begsample(i):endsample(i)) vdatend(1,endsample(i)) vdatend(2,endsample(i)) vdat(2,endsample(i):-1:begsample(i)) vdatbeg(2,begsample(i))];
      if isempty(color)
        color='none'; % this ensures that no edgelines will be plotted on top later on
      end
      h = patch(X, Y, [.6 .6 .6], 'FaceColor', facecolor, 'FaceAlpha', facealpha);
      set(h, 'linestyle', 'no');
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
      h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'markersize', markersize, 'markerfacecolor', markerfacecolor);
    elseif isequal(color, 'none')
      % do not plot the lines, this is useful in combination with highlightstyle=difference
      h = [];
    elseif ischar(color) && numel(color)==1
      % plot all lines with the same color
      h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'Color', color, 'markersize', markersize, 'markerfacecolor', markerfacecolor);
    elseif isnumeric(color) && isequal(size(color), [1 3])
      % plot all lines with the same RGB color
      h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'Color', color, 'markersize', markersize, 'markerfacecolor', markerfacecolor);
    elseif ischar(color) && numel(color)==nline
      % plot each line with its own color
      for i=1:size(vdat,1)
        h(i) = plot(hdat, vdat(i,:), style, 'LineWidth', linewidth, 'Color', color(i), 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      end
    elseif isnumeric(color) && size(color,1)==nline
      % the color is specified as Nx3 matrix with RGB values for each line
      for i=1:size(vdat,1)
        h(i) = plot(hdat, vdat(i,:), style, 'LineWidth', linewidth, 'Color', color(i,:), 'markersize', markersize, 'markerfacecolor', markerfacecolor);
      end
    elseif isnumeric(color) && size(color,1)==npos
      % the color is specified as Nx3 matrix with RGB values and varies over the length of the line
      for i=1:(size(vdat,2)-1)
        for j=1:size(vdat,1)
          h(i,j) = plot(hdat(i:i+1), vdat(j,i:i+1), style, 'LineWidth', linewidth, 'Color', mean(color([i i+1],:),1), 'markersize', markersize, 'markerfacecolor', markerfacecolor);
        end
      end
    else
      ft_warning('do not know how to plot the lines in the appropriate color');
      h = [];
    end
end % switch highlightstyle

if ~isempty(label)
  switch labelpos
    case 'upperleft'
      text(hpos - width/2, vpos + height/2, label, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
    case 'upperright'
      text(hpos + width/2, vpos + height/2, label, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
    case 'lowerleft'
      text(hpos - width/2, vpos - height/2, label, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
    case 'lowerright'
      text(hpos + width/2, vpos - height/2, label, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
    case 'center'
      text(hpos, vpos, label, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
  end
end

if box
  % this plots a box around the original hpos/vpos with appropriate width/height
  boxposition = zeros(1,4);
  boxposition(1) = hpos - width/2;
  boxposition(2) = hpos + width/2;
  boxposition(3) = vpos - height/2;
  boxposition(4) = vpos + height/2;
  ft_plot_box(boxposition);
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
    
    y_intercept = [0 0]; % If the y-axis crosses zero, the horizontal line should be at y = 0.
    
    % If the y-axis is all positive or negative, it should be as close to
    % zero as possible.
    if all(vlim > 0)
      y_intercept = repmat(min(vlim), 1, 2);
    elseif all(vlim < 0)
      y_intercept = repmat(max(vlim), 1, 2);
    end %if
    
    ft_plot_line(xrange, y_intercept, 'hpos', hpos, 'vpos', vpos, 'hlim', hlim, 'vlim', vlim, 'width', width, 'height', height);
  end
  if yaxis
    yrange = vlim;
    
    x_intercept = [0 0]; % If the x-axis crosses zero, the vertical line should be at x = 0.
    
    % If the x-axis is all positive or negative, it should be as close to
    % zero as possible.
    if all(hlim > 0)
      x_intercept = repmat(min(hlim), 1, 2);
    elseif all(vlim < 0)
      x_intercept = repmat(max(hlim), 1, 2);
    end %if
    
    ft_plot_line(x_intercept, yrange, 'hpos', hpos, 'vpos', vpos, 'hlim', hlim, 'vlim', vlim, 'width', width, 'height', height);
  end
end

set(h, 'tag', tag);

if ~isempty(parent)
  set(h, 'Parent', parent);
end

% the (optional) output is the handle
if nargout == 1
  varargout{1} = h;
end

if ~holdflag
  hold('off')
end

ft_warning(ws); % revert to original state
