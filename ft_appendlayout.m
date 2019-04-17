function [combined] = ft_appendlayout(cfg, varargin)

% FT_APPENDLAYOUT concatenates multiple layout descriptions that have been constructed
% separately.
%
% Use as
%   combined = ft_appendlayout(cfg, layout1, layout2, ...)
% where the input layouts result from FT_PREPARE_LAYOUT and the configuration
% should contain
%   cfg.direction = string, 'horizontal' or 'vertical' (default = 'horizontal')
%   cfg.align     = string, 'center', 'left', 'right', 'top' or 'bottom' (default = 'center')
%   cfg.distance  = number, distance between layouts (default is automatic)
%
% See also FT_PREPARE_LAYOUT, FT_LAYOUTPLOT, FT_APPENDSENS

% Undocumented options
%   cfg.xscale    = number, scaling to apply to input layouts along the horizontal direction (default = 1)
%   cfg.yscale    = number, scaling to apply to input layouts along the vertical direction (default = 1)

% Copyright (C) 2019, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% get the defaults
cfg.direction = ft_getopt(cfg, 'direction', 'horizontal');
cfg.align     = ft_getopt(cfg, 'align', 'center');
cfg.distance  = ft_getopt(cfg, 'distance', []);	% default will be determined further down
cfg.xscale    = ft_getopt(cfg, 'xscale', 1);
cfg.yscale    = ft_getopt(cfg, 'yscale', 1);

% check if the input data is valid for this function
% and ensure it is up to the latest standards
for i=1:length(varargin)
  assert(ft_datatype(varargin{i}, 'layout'), 'incorrect input, should be a layout structure');
end

% ensure these are cell-arrays
for i=1:numel(varargin)
  if isempty(varargin{i}.outline)
    varargin{i}.outline = {};
  end
  if isempty(varargin{i}.mask)
    varargin{i}.mask = {};
  end
end

% apply the scaling
for i=1:numel(varargin)
  varargin{i} = scalelayout(varargin{i}, cfg.xscale, cfg.yscale);
end

% remove the COMNT and SCALE
for i=1:numel(varargin)
  sel = ismember(varargin{i}.label, {'COMNT', 'SCALE'});
  varargin{i}.label  = varargin{i}.label(~sel);
  varargin{i}.pos    = varargin{i}.pos(~sel,:);
  varargin{i}.width  = varargin{i}.width(~sel);
  varargin{i}.height = varargin{i}.height(~sel);
end

% determine the center and size of each layout
for i=1:numel(varargin)
  xmin(i) = min(varargin{i}.pos(:,1) - varargin{i}.width/2);
  xmax(i) = max(varargin{i}.pos(:,1) + varargin{i}.width/2);
  ymin(i) = min(varargin{i}.pos(:,2) - varargin{i}.height/2);
  ymax(i) = max(varargin{i}.pos(:,2) + varargin{i}.height/2);
end
xrange  = xmax-xmin;
yrange  = ymax-ymin;
xcenter = (xmin+xmax)/2;
ycenter = (ymin+ymax)/2;
% extend these with a nan to facilitate the for-loop further down
xrange(end+1) = nan;
yrange(end+1) = nan;

% these are needed to automatically determine the distance
width = varargin{1}.width(:);
for i=1:numel(varargin)
  width = cat(1, width, varargin{i}.width(:));
end
height = varargin{1}.height(:);
for i=1:numel(varargin)
  width = cat(1, height, varargin{i}.height(:));
end

% set the first target point, it will be updated if we go along
xtarget = 0;
ytarget = 0;

% shift the anchor point of each of the layouts to its target position
switch cfg.direction
  case 'horizontal'
    if isempty(cfg.distance)
      % determine the distance between layouts
      distance = mean(width);
    else
      % use the specified distance
      distance = cfg.distance;
    end
    
    switch cfg.align
      case 'center'
        for i=1:numel(varargin)
          dx = xtarget - xcenter(i);
          dy = ytarget - ycenter(i);
          varargin{i} = shiftlayout(varargin{i}, dx, dy);
          xtarget = xtarget + xrange(i)/2 + xrange(i+1)/2 + distance;
        end % for varargin
        
      case 'top'
        for i=1:numel(varargin)
          dx = xtarget - xcenter(i);
          dy = ytarget - ymax(i);
          varargin{i} = shiftlayout(varargin{i}, dx, dy);
          xtarget = xtarget + xrange(i)/2 + xrange(i+1)/2 + distance;
        end % for varargin
        
      case 'bottom'
        for i=1:numel(varargin)
          dx = xtarget - xcenter(i);
          dy = ytarget - ymin(i);
          varargin{i} = shiftlayout(varargin{i}, dx, dy);
          xtarget = xtarget + xrange(i)/2 + xrange(i+1)/2 + distance;
        end % for varargin
        
      otherwise
        ft_error('invalid value for cfg.align');
    end % switch align
    
  case 'vertical'
    if isempty(cfg.distance)
      % determine the distance between layouts
      distance = mean(height);
    else
      % use the specified distance
      distance = cfg.distance;
    end
    
    switch cfg.align
      case 'center'
        for i=1:numel(varargin)
          dx = xtarget - xcenter(i);
          dy = ytarget - ycenter(i);
          varargin{i} = shiftlayout(varargin{i}, dx, dy);
          ytarget = ytarget - yrange(i)/2 - yrange(i+1)/2 - distance;
        end % for varargin
        
      case 'left'
        for i=1:numel(varargin)
          dx = xtarget - xmin(i);
          dy = ytarget - ycenter(i);
          varargin{i} = shiftlayout(varargin{i}, dx, dy);
          ytarget = ytarget - yrange(i)/2 - yrange(i+1)/2 - distance;
        end % for varargin
        
      case 'right'
        for i=1:numel(varargin)
          dx = xtarget - xmax(i);
          dy = ytarget - ycenter(i);
          varargin{i} = shiftlayout(varargin{i}, dx, dy);
          ytarget = ytarget - yrange(i)/2 - yrange(i+1)/2 - distance;
        end % for varargin
        
      otherwise
        ft_error('invalid value for cfg.align');
    end % switch align
    
  otherwise
    ft_error('invalid value for cfg.direction');
end % switch direction

% now that all layouts have been shifted, we can combine them
combined         = [];
combined.pos     = varargin{1}.pos;
combined.label   = varargin{1}.label(:);
combined.width   = varargin{1}.width(:);
combined.height  = varargin{1}.height(:);
combined.outline = varargin{1}.outline(:);
combined.mask    = varargin{1}.mask(:);
for i=2:numel(varargin)
  combined.pos     = cat(1, combined.pos,        varargin{i}.pos);
  combined.label   = cat(1, combined.label,      varargin{i}.label(:));
  combined.width   = cat(1, combined.width,      varargin{i}.width(:));
  combined.height  = cat(1, combined.height,     varargin{i}.height(:));
  combined.outline = cat(1, combined.outline(:), varargin{i}.outline(:));
  combined.mask    = cat(1, combined.mask(:),    varargin{i}.mask(:));
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance combined
ft_postamble history combined
ft_postamble savevar combined

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = shiftlayout(layout, dx, dy)
layout.pos(:,1) = layout.pos(:,1) + dx;
layout.pos(:,2) = layout.pos(:,2) + dy;
for i=1:numel(layout.mask)
  layout.mask{i}(:,1) = layout.mask{i}(:,1) + dx;
  layout.mask{i}(:,2) = layout.mask{i}(:,2) + dy;
end
for i=1:numel(layout.outline)
  layout.outline{i}(:,1) = layout.outline{i}(:,1) + dx;
  layout.outline{i}(:,2) = layout.outline{i}(:,2) + dy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = scalelayout(layout, sx, sy)
layout.pos(:,1) = layout.pos(:,1) * sx;
layout.pos(:,2) = layout.pos(:,2) * sy;
for i=1:numel(layout.mask)
  layout.mask{i}(:,1) = layout.mask{i}(:,1) * sx;
  layout.mask{i}(:,2) = layout.mask{i}(:,2) * sy;
end
for i=1:numel(layout.outline)
  layout.outline{i}(:,1) = layout.outline{i}(:,1) * sx;
  layout.outline{i}(:,2) = layout.outline{i}(:,2) * sy;
end
