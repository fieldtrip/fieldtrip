function ft_plot_lay(lay, varargin)

% FT_PLOT_LAY plots a two-dimensional layout
%
% Use as
%   ft_plot_lay(layout, ...)
% where the layout is a FieldTrip structure obtained from FT_PREPARE_LAYOUT.
%
% Additional options should be specified in key-value pairs and can be
%   'point'       = yes/no
%   'box'         = yes/no
%   'label'       = yes/no
%   'labelsize'   = number indicating font size (e.g. 6)
%   'labeloffset' = offset of label from point (suggestion is 0.005)
%   'labelrotate' = scalar, vector with rotation angle (in degrees) per label (default = 0)
%   'labelalignh' = string, or cell-array specifying the horizontal alignment of the text (default = 'left')
%   'labelalignv' = string, or cell-array specifying the vertical alignment of the text (default = 'middle')
%   'mask'        = yes/no
%   'outline'     = yes/no
%   'verbose'     = yes/no
%   'pointsymbol' = string with symbol (e.g. 'o') - all three point options need to be used together
%   'pointcolor'  = string with color (e.g. 'k')
%   'pointsize'   = number indicating size (e.g. 8)
%
% It is possible to plot the object in a local pseudo-axis (c.f. subplot), which is specfied as follows
%   'hpos'        = horizontal position of the lower left corner of the local axes
%   'vpos'        = vertical position of the lower left corner of the local axes
%   'width'       = width of the local axes
%   'height'      = height of the local axes

% Copyright (C) 2009, Robert Oostenveld
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
hpos         = ft_getopt(varargin, 'hpos',         0);
vpos         = ft_getopt(varargin, 'vpos',         0);
width        = ft_getopt(varargin, 'width',          []);
height       = ft_getopt(varargin, 'height',         []);
point        = ft_getopt(varargin, 'point',        true);
box          = ft_getopt(varargin, 'box',          true);
label        = ft_getopt(varargin, 'label',        true);
labelsize    = ft_getopt(varargin, 'labelsize',    10);
labelfont    = ft_getopt(varargin, 'labelfont',    'helvetica');
labeloffset  = ft_getopt(varargin, 'labeloffset',  0);
mask         = ft_getopt(varargin, 'mask',         true);
outline      = ft_getopt(varargin, 'outline',      true);
verbose      = ft_getopt(varargin, 'verbose',      false);
pointsymbol  = ft_getopt(varargin, 'pointsymbol');
pointcolor   = ft_getopt(varargin, 'pointcolor');
pointsize    = ft_getopt(varargin, 'pointsize');
interpreter  = ft_getopt(varargin, 'interpreter', 'tex');

% some stuff related to some refined label plotting
labelrotate   = ft_getopt(varargin, 'labelrotate',  0);
labelalignh   = ft_getopt(varargin, 'labelalignh',  'left');
labelalignv   = ft_getopt(varargin, 'labelalignv',  'middle');
labelcolor    = ft_getopt(varargin, 'labelcolor', 'k');


% convert between true/false/yes/no etc. statements
point   = istrue(point);
box     = istrue(box);
label   = istrue(label);
mask    = istrue(mask);
outline = istrue(outline);
verbose = istrue(verbose);

if ~(point || box || label || mask || outline)
  % there is nothing to be plotted
  return;
end

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

% layout units can be arbitrary (e.g. pixels for .mat files)
% so we need to compute the right scaling factor and offset
% create a matrix with all coordinates
% from positions, mask, and outline
allCoords = lay.pos;
if isfield(lay, 'mask') && ~isempty(lay.mask)
  for k = 1:numel(lay.mask)
    allCoords = [allCoords; lay.mask{k}];
  end
end
if isfield(lay, 'outline') &&~isempty(lay.outline)
  for k = 1:numel(lay.outline)
    allCoords = [allCoords; lay.outline{k}];
  end
end

naturalWidth = (max(allCoords(:,1))-min(allCoords(:,1)));
naturalHeight = (max(allCoords(:,2))-min(allCoords(:,2)));

if isempty(width) && isempty(height)
  xScaling = 1;
  yScaling = 1;
elseif isempty(width) && ~isempty(height)
  % height specified, auto-compute width while maintaining aspect ratio
  yScaling = height/naturalHeight;
  xScaling = yScaling;
elseif ~isempty(width) && isempty(height)
  % width specified, auto-compute height while maintaining aspect ratio
  xScaling = width/naturalWidth;
  yScaling = xScaling;
else
  % both width and height specified
  xScaling = width/naturalWidth;
  yScaling = height/naturalHeight;
end

X      = lay.pos(:,1)*xScaling + hpos;
Y      = lay.pos(:,2)*yScaling + vpos;
Width  = lay.width*xScaling;
Height = lay.height*yScaling;
Lbl    = lay.label;

if point
  if ~isempty(pointsymbol) && ~isempty(pointcolor) && ~isempty(pointsize) % if they're all non-empty, don't use the default
    plot(X, Y, 'marker',pointsymbol,'color',pointcolor,'markersize',pointsize,'linestyle','none');
  else
    plot(X, Y, 'marker','.','color','b','linestyle','none');
    plot(X, Y, 'marker','o','color','y','linestyle','none');
  end
end

if label
  % the MATLAB text function fails if the position for the string is specified in single precision
  X = double(X);
  Y = double(Y);
  
  % check whether fancy label plotting is needed, this requires a for loop,
  % otherwise print text in a single shot
  if numel(labelrotate)==1
    text(X+labeloffset, Y+(labeloffset*1.5), Lbl ,'fontsize',labelsize,'fontname',labelfont,'interpreter',interpreter,'horizontalalignment',labelalignh,'verticalalignment',labelalignv,'color',labelcolor);
  else
    n = numel(Lbl);
    if ~iscell(labelalignh)
      labelalignh = repmat({labelalignh},[n 1]);
    end
    if ~iscell(labelalignv)
      labelalignv = repmat({labelalignv},[n 1]);
    end
    if numel(Lbl)~=numel(labelrotate)||numel(Lbl)~=numel(labelalignh)||numel(Lbl)~=numel(labelalignv)
      eror('there is something wrong with the input arguments');
    end
    for k = 1:numel(Lbl)
      text(X(k)+labeloffset, Y(k)+(labeloffset*1.5), Lbl{k}, 'fontsize', labelsize, 'fontname', labelfont, 'interpreter', interpreter, 'horizontalalignment', labelalignh{k}, 'verticalalignment', labelalignv{k}, 'rotation', labelrotate(k),'color',labelcolor);
    end
  end
end

if box
  line([X-Width/2 X+Width/2 X+Width/2 X-Width/2 X-Width/2]',[Y-Height/2 Y-Height/2 Y+Height/2 Y+Height/2 Y-Height/2]', 'color', [0 0 0]);
end

if outline && isfield(lay, 'outline')
  if verbose
    fprintf('solid lines indicate the outline, e.g. head shape or sulci\n');
  end
  for i=1:length(lay.outline)
    if ~isempty(lay.outline{i})
      X = lay.outline{i}(:,1)*xScaling + hpos;
      Y = lay.outline{i}(:,2)*yScaling + vpos;
      h = line(X, Y);
      set(h, 'color', 'k');
      set(h, 'linewidth', 2);
    end
  end
end

if mask && isfield(lay, 'mask')
  if verbose
    fprintf('dashed lines indicate the mask for topograpic interpolation\n');
  end
  for i=1:length(lay.mask)
    if ~isempty(lay.mask{i})
      X = lay.mask{i}(:,1)*xScaling + hpos;
      Y = lay.mask{i}(:,2)*yScaling + vpos;
      % the polygon representing the mask should be closed
      X(end+1) = X(1);
      Y(end+1) = Y(1);
      h = line(X, Y);
      set(h, 'color', 'k');
      set(h, 'linewidth', 1.5);
      set(h, 'linestyle', ':');
    end
  end
end

axis auto
axis equal
axis off

if ~holdflag
  hold off
end

warning(ws); %revert to original state
