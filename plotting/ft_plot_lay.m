function ft_plot_lay(lay, varargin)

% FT_PLOT_LAY plots a two-dimensional layout
%
% Use as
%   ft_plot_lay(layout, ...)
% where the layout is a FieldTrip structure obtained from FT_PREPARE_LAYOUT.
%
% Additional options should be specified in key-value pairs and can be
%   'point'         = yes/no
%   'box'           = yes/no
%   'label'         = yes/no
%   'labelsize'     = number indicating font size (e.g. 6)
%   'labeloffset'   = offset of label from point (suggestion is 0.005)
%   'mask'          = yes/no
%   'outline'       = yes/no
%   'verbose'       = yes/no
%   'pointsymbol'   = string with symbol (e.g. 'o') - all three point options need to be used together
%   'pointcolor'    = string with color (e.g. 'k')
%   'pointsize'     = number indicating size (e.g. 8)

% Copyright (C) 2009, Robert Oostenveld
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

warning('on', 'MATLAB:divideByZero');

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'point', 'box', 'label','labelsize','labeloffset', 'mask', 'outline', 'verbose','pointsymbol','pointcolor','pointsize'});
hpos         = keyval('hpos',        varargin); if isempty(hpos),       hpos = 0;            end
vpos         = keyval('vpos',        varargin); if isempty(vpos),       vpos = 0;            end
point        = keyval('point',       varargin); if isempty(point),      point = true;        end
box          = keyval('box',         varargin); if isempty(box),        box = true;          end
label        = keyval('label',       varargin); if isempty(label),      label = true;        end
labelsize    = keyval('labelsize',   varargin); if isempty(labelsize),  labelsize = 10;      end
labeloffset  = keyval('labeloffset', varargin); if isempty(labeloffset),labeloffset = 0;     end
mask         = keyval('mask',        varargin); if isempty(mask),       mask = true;         end
outline      = keyval('outline',     varargin); if isempty(outline),    outline = true;      end
verbose      = keyval('verbose',     varargin); if isempty(verbose),    verbose = false;     end
pointsymbol  = keyval('pointsymbol', varargin);
pointcolor   = keyval('pointcolor',  varargin);
pointsize    = keyval('pointsize',   varargin);

% convert between true/false/yes/no etc. statements
point   = istrue(point);
box     = istrue(box);
label   = istrue(label);
mask    = istrue(mask);
outline = istrue(outline);
verbose = istrue(verbose);


% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

X      = lay.pos(:,1) + hpos;
Y      = lay.pos(:,2) + vpos;
Width  = lay.width;
Height = lay.height;
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
  text(X+labeloffset, Y+(labeloffset*1.5), Lbl,'fontsize',labelsize);
end

if box
  line([X-Width/2 X+Width/2 X+Width/2 X-Width/2 X-Width/2]',[Y-Height/2 Y-Height/2 Y+Height/2 Y+Height/2 Y-Height/2]');
end

if outline && isfield(lay, 'outline')
  if verbose  
    fprintf('solid lines indicate the outline, e.g. head shape or sulci\n');
  end
  for i=1:length(lay.outline)
    if ~isempty(lay.outline{i})
      X = lay.outline{i}(:,1) + hpos;
      Y = lay.outline{i}(:,2) + vpos;
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
      X = lay.mask{i}(:,1) + hpos;
      Y = lay.mask{i}(:,2) + vpos;
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
