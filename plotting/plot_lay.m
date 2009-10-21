function plot_lay(lay, varargin)

% PLOT_LAY plots a two-dimensional layout
%
% Use as
%   plot_lay(layout, ...)
% where the layout is a FieldTrip structure obtained from PREPARE_LAYOUT.
%
% Additional options should be specified in key-value pairs and can be
%   'point'    = yes/no
%   'box'      = yes/no
%   'label'    = yes/no
%   'mask'     = yes/no
%   'outline'  = yes/no

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: plot_lay.m,v $
% Revision 1.3  2009/06/02 15:46:37  giopia
% added vpos and hpos for position and verbose to control fprintf
%
% Revision 1.2  2009/05/12 12:34:27  roboos
% updated docuemntation
%
% Revision 1.1  2009/05/12 12:30:19  roboos
% first version, used by Chris for interactive cohrefchan picking
%

point   = keyval('point',   varargin{:}); if isempty(point),    point = true;   end
box     = keyval('box',     varargin{:}); if isempty(box),      box = true;     end
label   = keyval('label',   varargin{:}); if isempty(label),    label = true;   end
mask    = keyval('mask',    varargin{:}); if isempty(mask),     mask = true;    end
outline = keyval('outline', varargin{:}); if isempty(outline),  outline = true; end
hpos    = keyval('hpos',    varargin{:}); if isempty(hpos),     hpos = 0;       end
vpos    = keyval('vpos',    varargin{:}); if isempty(vpos),     vpos = 0;       end
verbose = keyval('verbose', varargin{:}); if isempty(verbose),  verbose = true; end

% everything is added to the current figure
holdflag = ishold;
hold on

X      = lay.pos(:,1) + hpos;
Y      = lay.pos(:,2) + vpos;
Width  = lay.width;
Height = lay.height;
Lbl    = lay.label;

if point
  plot(X, Y, 'b.');
  plot(X, Y, 'yo');
end

if label
  text(X, Y, Lbl);
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
