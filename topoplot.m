function [handle] = topoplot(varargin)

% TOPOPLOT plots a topographic map of an EEG or MEG field as a 2-D
% circular view (looking down at the top of the head) using interpolation
% on a fine cartesian grid.
%
% This function is called by topoplotER or topoplotTFR
%
% You can also call this function directly as follows:
%         topoplot(cfg, datavector)
%         topoplot(cfg, X, Y, datavector)
%         topoplot(cfg, X, Y, datavector, Labels)
%         topoplot(datavector,'Key1','Value1','Key2','Value2',...)
%
% Inputs can be either:
%     datavector  = vector of values to be plotted as color
%     cfg         = configuration structure containing the (optional)
%     parameters
%     X           = x-coordinates for channels in datavector
%     Y           = y-coordinates for channels in datavector
%     Labels      = labels for channels in datavector
% or the inputs can be key-value pairs containing the (optional) parameters.
% Every cfg field can be specified using the fieldname as a key in the
% key-value pairs.
%
% if X, Y and Labels are given, cfg.layout is NOT used. If X, Y, and Labels
% are not given, cfg.layout must be given and it is assumed that the
% channels in the datavector exactly mach the channels in the layout.
%
% The layout defines how the channels will be arranged in the 2-D plane.
% You can specify the layout in a variety of ways:
%  - you can give the name of an ascii layout file with extension *.lay
%  - you can give the name of an electrode file
%  - you can give an electrode definition,  i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% If you do not specify any of these, and if the data structure contains an
% electrode or gradiometer structure, that will be used for creating a
% layout.
%
% Optional Parameters and Values
%
% cfg.colormap        = any sized colormap, see COLORMAP
% cfg.colorbar        = 'yes'
%                       'no' (default)
%                       'North'              inside plot box near top
%                       'South'              inside bottom
%                       'East'               inside right
%                       'West'               inside left
%                       'NorthOutside'       outside plot box near top
%                       'SouthOutside'       outside bottom
%                       'EastOutside'        outside right
%                       'WestOutside'        outside left
% cfg.interplimits    = limits for interpolation (default = 'head')
%                       'electrodes' to furthest electrode
%                       'head' to edge of head
% cfg.gridscale       = scaling grid size (default = 67)
%                       determines resolution of figure
% cfg.maplimits       = 'absmax' +/- the absolute-max (default = 'absmax')
%                       'maxmin' scale to data range
%                       [clim1, clim2] user-defined lo/hi
% cfg.style           = topoplot style (default = 'both')
%                       'straight' colormap only
%                       'contour' contour lines only
%                       'both' (default) both colormap and contour lines
%                       'fill' constant color between lines
%                       'blank' just head and electrodes
% cfg.contournum      = number of contour lines (default = 6), see CONTOUR
% cfg.shading         = 'flat' 'interp' (default = 'flat')
% cfg.interpolation   = 'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
% cfg.headcolor       = Color of head cartoon (default = [0,0,0])
% cfg.hlinewidth      = number, Linewidth of the drawn head, nose and ears (default = 2)
% cfg.contcolor       = Contourline color (default = [0 0 0])
% cfg.electrodes      = 'on','off','labels','numbers','highlights' or 'dotnum' (default = 'on')
% cfg.emarker         = Marker symbol (default = 'o')
% cfg.ecolor          = Marker color (default = [0 0 0] (black))
% cfg.emarkersize     = Marker size (default = 2)
% cfg.efontsize       = Font size of electrode labels/numbers (default = 8 pt)
%                       when cfg.electrodes = 'numbers' or 'labels'
% cfg.comment         =  string of text
% cfg.commentpos      = position of comment (default = 'leftbottom')
%                       'lefttop' 'leftbottom' 'middletop' 'middlebottom' 'righttop' 'rightbottom'
%                       or [x y] coordinates
%                       or 'title' to place comment as title
% cfg.fontsize        = Font size of comment (default = 8 pt)
% cfg.highlight       = 'off' or the channel numbers you want to highlight (default = 'off').
%                       These numbers should correspond with the channels in the data, not in
%                       the layout file.
% cfg.hlmarker        = Highlight marker symbol (default = 'o')
% cfg.hlcolor         = Highlight marker color (default = [0 0 0] (black))
% cfg.hlmarkersize    = Highlight marker size (default = 6)
% cfg.hllinewidth     = Highlight marker linewidth (default = 3)
% cfg.outline         = 'scalp' or 'ECog' (default = 'scalp')
%
% Note: topoplot() only works when map limits are >= the max and min
%                             interpolated data values.

% Undocumented local options:
% cfg.grid
% cfg.maxchans
% cfg.showlabels
% cfg.zlim
% cfg.mask for opacity masking, e.g. with statistical significance

% Copyright (C) 1996, Andy Spydell, Colin Humphries & Arnaud Delorme, CNL / Salk Institute
% Copyright (C) 2004-2009, F.C. Donders Centre, New implementation by Geerten Kramer, based on versions of Ole Jensen and Jan-Mathijs Schoffelen
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

% Try to detect EEGLAB-style input and give an informative error
% message. The EEGLAB documentation describes the usage as
%        >>  topoplot(datavector, EEG.chanlocs);   % plot a map using an EEG chanlocs structure
%        >>  topoplot(datavector, 'my_chan.locs'); % read a channel locations file and plot a map
%        >>  topoplot('example');                  % give an example of an electrode location file
%        >>  [h grid_or_val plotrad_or_grid, xmesh, ymesh]= ...
%                           topoplot(datavector, chan_locs, 'Input1','Value1', ...);

if nargin==2 && isvector(varargin{1}) && isstruct(varargin{2}) && isfield(varargin{2}, 'labels')
  eeglab = 1;
elseif nargin==2 && isvector(varargin{1}) && ischar(varargin{2})
  eeglab = 1;
elseif nargin==1 && isequal(varargin{1}, 'example')
  eeglab = 1;
elseif nargin>2 && isvector(varargin{1}) && mod(nargin,2)==0 && isstruct(varargin{2}) && isfield(varargin{2}, 'labels')
  eeglab = 1;
else
  eeglab = 0;
end

if eeglab
  % the input resembles the input of the EEGLAB topoplot function
  error('Unrecognized input, please look at "help topoplot", "which topoplot" and "path". The input looks as if you expect the EEGLAB version of the topoplot function. Your path settings may be incorrect and the FieldTrip and EEGLAB version of the "topoplot" function may be confused.')
end

% deal with the different types of input syntax that this function can get
if mod(nargin,2) && isnumeric(varargin{1}) && ischar(varargin{2})
  % topoplot(data, key, val, ...)
  cfg  = keyval2cfg(varargin(2:end));
  data = varargin{1};
  OldStyleSyntax=0;
elseif nargin==2
  % topoplot(cfg,data)
  OldStyleSyntax=0;
  cfg  = varargin{1};
  data = varargin{2};
  err  = 0;
  if ~isempty(cfg),
    err  = err + double(~isstruct(cfg) && ~strcmp(class(cfg), 'config'));
  end
  err  = err + ~isnumeric(data);
  if err
    errmsg=['\n'];
    errmsg=[errmsg,'When two input arguments are supplied, the following syntax should be used:\n'];
    errmsg=[errmsg,'topoplot(cfg,datavector);\n'];
    error(sprintf(errmsg));
  end;
elseif nargin==4
  % topoplot(cfg,X,Y,data)
  OldStyleSyntax=1;
  cfg   = varargin{1};
  chanX = varargin{2};
  chanY = varargin{3};
  data  = varargin{4};
  err   = 0;
  if ~isempty(cfg),
    err  = err + double(~isstruct(cfg) && ~strcmp(class(cfg), 'config'));
  end
  err  = err + ~isnumeric(data);
  err  = err + ~isnumeric(chanX);
  err  = err + ~isnumeric(chanY);
  if err
    errmsg=['\n'];
    errmsg=[errmsg,'When four input arguments are supplied, the following syntax should be used:\n'];
    errmsg=[errmsg,'topoplot(cfg,X,Y,datavector);\n'];
    error(sprintf(errmsg));
  end;
elseif nargin==5
  % topoplot(cfg,X,Y,data,labels)
  OldStyleSyntax=1;
  cfg    = varargin{1};
  chanX  = varargin{2};
  chanY  = varargin{3};
  data   = varargin{4};
  chanLabels = varargin{5};
  err    = 0;
  if ~isempty(cfg),
    err  = err + double(~isstruct(cfg) && ~strcmp(class(cfg), 'config'));
  end
  err    = err + ~isnumeric(data);
  err    = err + ~isnumeric(chanX);
  err    = err + ~isnumeric(chanY);
  err    = err + ~iscell(chanLabels);
  err    = err + numel(chanLabels)~=numel(chanX);
  err    = err + numel(chanLabels)~=numel(chanY);
  if err
    errmsg=['\n'];
    errmsg=[errmsg,'When five input arguments are supplied, the following syntax should be used:\n'];
    errmsg=[errmsg,'topoplot(cfg,X,Y,datavector,Labels);\n'];
    error(sprintf(errmsg));
  end;
else
  error('unrecognized input, please look at the help of this function')
end

% set the defaults
if ~isfield(cfg, 'maxchans')      cfg.maxchans = 256;       end;
if ~isfield(cfg, 'maplimits')     cfg.maplimits = 'absmax'; end; % absmax, maxmin, [values]
if ~isfield(cfg, 'interplimits')  cfg.interplimits ='head'; end; % head, electrodes
if ~isfield(cfg, 'gridscale')     cfg.gridscale = 67;       end; % 67 in original
if ~isfield(cfg, 'contournum')    cfg.contournum = 6;       end;
if ~isfield(cfg, 'colorbar')      cfg.colorbar = 'no';      end;
if ~isfield(cfg, 'style')         cfg.style = 'both';       end; % both,straight,fill,contour,blank
if ~isfield(cfg, 'headcolor')     cfg.headcolor = [0 0 0];  end;
if ~isfield(cfg, 'contcolor')     cfg.contcolor = 'k';      end;
if ~isfield(cfg, 'hlinewidth')    cfg.hlinewidth = 2;       end;
if ~isfield(cfg, 'shading')       cfg.shading = 'flat';     end; % flat or interp
if ~isfield(cfg, 'interpolation') cfg.interpolation = 'v4'; end;
if ~isfield(cfg, 'fontsize'),     cfg.fontsize = 8;         end;
if ~isfield(cfg, 'commentpos'),   cfg.commentpos = 'leftbottom';    end;
if ~isfield(cfg, 'mask'),         cfg.mask = [];            end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to update the existing figure, this is to speed up realtime plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(cfg, 'update') && strcmp(cfg.update, 'yes')
  update = guidata(gcf);
  if ~isempty(update)
    Xi = update.Xi;
    Yi = update.Yi;
    x  = update.x;
    y  = update.y;
    xi = update.xi;
    yi = update.yi;
    % Interpolate the topographic data
    Zi = griddata(x', y, data, xi', yi, cfg.interpolation);
    % keep the same NaNs
    Zi(isnan(update.Zi)) = NaN;
    deltax = xi(2)-xi(1); % length of grid entry
    deltay = yi(2)-yi(1); % length of grid entry
    surface(Xi-deltax/2,Yi-deltay/2,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',cfg.shading);
    return
  end
end

if ~isfield(cfg,'layout')
  if ~OldStyleSyntax
    error('Specify at least the field or key "layout".');
  end;
end;

if ~ischar(cfg.contcolor)     cfg.contcolor = 'k'; warning('cfg.contcolor must be string, put to ''k''');   end;

if ~isfield(cfg,'electrodes')   cfg.electrodes = 'on'; end; % on,off,label,numbers or highlights
if ~isfield(cfg,'showlabels') % for compatibility with OLDSTYLE
  cfg.showlabels = '';
else
  cfg.electrodes = '';
end;

if ~isfield(cfg,'emarker')      cfg.emarker = 'o';     end;
if ~isfield(cfg,'ecolor')       cfg.ecolor = [0 0 0];  end;
if ~isfield(cfg,'emarkersize')  cfg.emarkersize = 2;   end;
if ~isfield(cfg,'efontsize')    cfg.efontsize = get(0,'DefaultAxesFontSize');end;

if ~isfield(cfg,'highlight')    cfg.highlight = 'off'; end; % 'off' or the electrodenumbers.
if ~isfield(cfg,'hlmarker')     cfg.hlmarker = 'o';    end;
if ~isfield(cfg,'hlcolor')      cfg.hlcolor = [0 0 0]; end;
if ~isfield(cfg,'hlmarkersize') cfg.hlmarkersize = 6;  end;
if ~isfield(cfg,'hllinewidth')  cfg.hllinewidth = 3;   end;
if ~isfield(cfg,'hlfacecolor')  cfg.hlfacecolor = cfg.hlcolor; end;

if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('topoplot(): Colormap must be a n x 3 matrix'); end
  colormap(cfg.colormap);
end;

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'renamed',     {'grid_scale',  'gridscale'});
cfg = checkconfig(cfg, 'renamed',     {'interpolate', 'interpolation'});
cfg = checkconfig(cfg, 'renamed',     {'numcontour',  'contournum'});
cfg = checkconfig(cfg, 'renamed',     {'electrod',    'electrodes'});
cfg = checkconfig(cfg, 'renamed',     {'hcolor',      'headcolor'});
cfg = checkconfig(cfg, 'renamed',     {'electcolor',  'ecolor'});
cfg = checkconfig(cfg, 'renamed',     {'emsize',      'emarkersize'});
cfg = checkconfig(cfg, 'renamed',     {'efsize',      'efontsize'});
cfg = checkconfig(cfg, 'renamed',     {'headlimits',  'interplimits'});

if isfield(cfg,'interplimits')
  if ~ischar(cfg.interplimits), error('topoplot(): interplimits value must be a string'); end
  cfg.interplimits = lower(cfg.interplimits);
  if ~strcmp(cfg.interplimits,'electrodes') && ~strcmp(cfg.interplimits,'head') && ~strcmp(cfg.interplimits,'headleft'),
    error('topoplot(): Incorrect value for interplimits');
  end
end;

if isfield(cfg,'shading')
  cfg.shading = lower(cfg.shading);
  if ~any(strcmp(cfg.shading,{'flat','interp'})), error('Invalid Shading Parameter'); end
end

% for compatibility with topoplotXXX functions
if isfield(cfg,'zlim')
  cfg.maplimits = cfg.zlim;
  cfg           = rmfield(cfg,'zlim');
end;

[numChan,numTime] = size(data);
if numChan && numTime>1
  error('topoplot(): data should be a column vector\n');
end

if ~OldStyleSyntax
  % create layout including channel positions, labels and anatomical mask and outline
  cfg.layout  = prepare_layout(cfg);
  chanX       = cfg.layout.pos(:,1);
  chanY       = cfg.layout.pos(:,2);
  chanLabels  = cfg.layout.label(:);
else
  % create layout including channel positions, labels and anatomical mask and outline
  cfg.layout.pos   = [chanX chanY];
  cfg.layout.label = chanLabels;
  cfg.layout       = prepare_layout(cfg);
end
clear chanX chanY

% The whole figure will be created with the x-axis along the horizontal
% dimension of the figure and the y-axis along the vertical dimension. This
% means that for data coming from a system with a head-coordinate system in
% which the x-axis points to the nose (e.g. CTF), the layout should be 90
% degrees rotated relative to the head coordinates.
x = cfg.layout.pos(:,1);
y = cfg.layout.pos(:,2);

if isfield(cfg, 'outline')
  error('the option cfg.outline is not supported any more, please contact Robert for a detailled explanation');
end

if exist('chanLabels', 'var'),
  ind_SCALE = strmatch('SCALE', chanLabels);
  if length(ind_SCALE)==1
    % remember the position of  the scale
    X_SCALE               = cfg.layout.pos(ind_SCALE, 1);
    Y_SCALE               = cfg.layout.pos(ind_SCALE, 2);
    x(ind_SCALE) = [];
    y(ind_SCALE) = [];
    chanLabels(ind_SCALE) = [];
  end
  ind_COMNT = strmatch('COMNT', chanLabels);
  if length(ind_COMNT)==1
    % remember the position of the comment
    X_COMNT               = cfg.layout.pos(ind_COMNT, 1);
    Y_COMNT               = cfg.layout.pos(ind_COMNT, 1);
    x(ind_COMNT) = [];
    y(ind_COMNT) = [];
    chanLabels(ind_COMNT) = [];
  end
end

% Set coordinates for comment
if strcmp(cfg.commentpos,'lefttop')
  x_COMNT = -0.7;
  y_COMNT =  0.6;
  HorAlign = 'left';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos,'leftbottom')
  x_COMNT = -0.6;
  y_COMNT = -0.6;
  HorAlign = 'left';
  VerAlign = 'bottom';
elseif strcmp(cfg.commentpos,'middletop')
  x_COMNT =  0;
  y_COMNT =  0.75;
  HorAlign = 'center';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos,'middlebottom')
  x_COMNT =  0;
  y_COMNT = -0.7;
  HorAlign = 'center';
  VerAlign = 'bottom';
elseif strcmp(cfg.commentpos,'righttop')
  x_COMNT =  0.65;
  y_COMNT =  0.6;
  HorAlign = 'right';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos,'rightbottom')
  x_COMNT =  0.6;
  y_COMNT = -0.6;
  HorAlign = 'right';
  VerAlign = 'bottom';
elseif isnumeric(cfg.commentpos)
  x_COMNT = cfg.commentpos(1);
  y_COMNT = cfg.commentpos(2);
  HorAlign = 'left';
  VerAlign = 'middle';
  x_COMNT = 0.9*((x_COMNT-min(x))/(max(x)-min(x))-0.5);
  y_COMNT = 0.9*((y_COMNT-min(y))/(max(y)-min(y))-0.5);
end

gca;
cla;
hold on

if ~strcmp(cfg.style,'blank')
  % find limits for interpolation:
  if strcmp(cfg.interplimits,'head')
    xmin = +inf;
    xmax = -inf;
    ymin = +inf;
    ymax = -inf;
    for i=1:length(cfg.layout.mask)
      xmin = min([xmin; cfg.layout.mask{i}(:,1)]);
      xmax = max([xmax; cfg.layout.mask{i}(:,1)]);
      ymin = min([ymin; cfg.layout.mask{i}(:,2)]);
      ymax = max([ymax; cfg.layout.mask{i}(:,2)]);
    end
  elseif strcmp(cfg.interplimits,'electrodes')
    xmin = min(cfg.layout.pos(:,1));
    xmax = max(cfg.layout.pos(:,1));
    ymin = min(cfg.layout.pos(:,2));
    ymax = max(cfg.layout.pos(:,2));
  elseif strcmp(cfg.interplimits,'headleft')
    xmin = +inf;
    xmax = -inf;
    ymin = +inf;
    ymax = -inf;
    for i=1:length(cfg.layout.mask)
      xmin = min([xmin; cfg.layout.mask{i}(:,1)]);
      xmax = 0.02;
      ymin = min([ymin; cfg.layout.mask{i}(:,2)]);
      ymax = max([ymax; cfg.layout.mask{i}(:,2)]);
    end
  else
    xmin = min(cfg.layout.pos(:,1));
    xmax = max(cfg.layout.pos(:,1));
    ymin = min(cfg.layout.pos(:,2));
    ymax = max(cfg.layout.pos(:,2));
  end

  xi         = linspace(xmin,xmax,cfg.gridscale);   % x-axis for interpolation (row vector)
  yi         = linspace(ymin,ymax,cfg.gridscale);   % y-axis for interpolation (row vector)
  [Xi,Yi,Zi] = griddata(x', y, data, xi', yi, cfg.interpolation); % Interpolate the topographic data

  % calculate colormap limits
  m = size(colormap,1);
  if ischar(cfg.maplimits)
    if strcmp(cfg.maplimits,'absmax')
      amin = -max(max(abs(Zi)));
      amax = max(max(abs(Zi)));
    elseif strcmp(cfg.maplimits,'maxmin')
      amin = min(min(Zi));
      amax = max(max(Zi));
    end
  else
    amin = cfg.maplimits(1);
    amax = cfg.maplimits(2);
  end
  % FIXME, according to Ingrid and Robert (19 Oct 2009), these deltas probably should be 0
  deltax = xi(2)-xi(1); % length of grid entry
  deltay = yi(2)-yi(1); % length of grid entry

  if isfield(cfg.layout, 'mask') && ~isempty(cfg.layout.mask)
    % apply anatomical mask to the data, i.e. that determines that the interpolated data outside the circle is not displayed
    maskA = false(size(Zi));
    for i=1:length(cfg.layout.mask)
      cfg.layout.mask{i}(end+1,:) = cfg.layout.mask{i}(1,:); % force them to be closed
      maskA(inside_contour([Xi(:) Yi(:)], cfg.layout.mask{i})) = true;
    end
    Zi(~maskA) = NaN;
  end

  if ~isempty(cfg.mask),
    % this mask is based on some statistical feature of the data itself, e.g. significance and is not related to the anatomical mask
    [maskX,maskY,maskZ] = griddata(x', y, double(cfg.mask), xi', yi, cfg.interpolation);
    % mask should be scaled between 0 and 1, clip the values that ly outside that range
    maskZ(isnan(maskZ)) = 0;
    maskZ(isinf(maskZ)) = 0;
    maskZ(maskZ<0) = 0;
    maskZ(maskZ>1) = 1;
  end

  % Draw topoplot on head
  if strcmp(cfg.style,'contour')
    contour(Xi,Yi,Zi,cfg.contournum,cfg.contcolor);
  elseif strcmp(cfg.style,'both')
    % first draw the surface, then the contour, to ensure that after exporting the contour lines are "on top"
    h = surface(Xi-deltax/2,Yi-deltay/2,zeros(size(Zi)),Zi,'EdgeColor','none', 'FaceColor',cfg.shading);
    if exist('maskZ','var'),
      set(h, 'AlphaData', maskZ);
      alim([0 1]);
      set(h, 'FaceAlpha', 'interp');
    end
    contour(Xi,Yi,Zi,cfg.contournum,cfg.contcolor);
  elseif strcmp(cfg.style,'straight')
    h = surface(Xi-deltax/2,Yi-deltay/2,zeros(size(Zi)),Zi,'EdgeColor','none', 'FaceColor',cfg.shading);
    if exist('maskZ','var'),
      set(h, 'AlphaData', maskZ);
      alim([0 1]);
      set(h, 'FaceAlpha', 'interp');
    end
  elseif strcmp(cfg.style,'fill')
    contourf(Xi,Yi,Zi,cfg.contournum,cfg.contcolor);
  else
    error('Invalid style')
  end
  caxis([amin amax]); % set coloraxis
end

% Plot electrodes:
if strcmp(cfg.electrodes,'on') || strcmp(cfg.showlabels,'markers')
  if ischar(cfg.highlight)
    plot_vector(x,y,'style',cfg.emarker,'Color',cfg.ecolor,'markersize',cfg.emarkersize);
  elseif isnumeric(cfg.highlight)
    normal = setdiff(1:length(x), cfg.highlight);
    plot_vector(x(normal), y(normal), 'style', cfg.emarker, 'Color', cfg.ecolor, 'markersize', cfg.emarkersize);
    plot_vector(x(cfg.highlight), y(cfg.highlight), 'style', cfg.hlmarker, 'Color', cfg.hlcolor, 'markersize', cfg.hlmarkersize, ...
      'linewidth',  cfg.hllinewidth, 'markerfacecolor', cfg.hlfacecolor);    
  elseif iscell(cfg.highlight)
    plot_vector(x,y,'style',cfg.emarker,'Color',cfg.ecolor,'markersize',cfg.emarkersize);
    for iCell = 1:length(cfg.highlight)
      plot_vector(x(cfg.highlight{iCell}), y(cfg.highlight{iCell}), 'style', cfg.hlmarker{iCell}, 'Color', ...
        cfg.hlcolor{iCell}, 'markersize', cfg.hlmarkersize{iCell}, ...
      'linewidth',  cfg.hllinewidth{iCell}, 'markerfacecolor', cfg.hlfacecolor{iCell});  
    end
  else
    error('Unknown highlight type');
  end;
elseif any(strcmp(cfg.electrodes,{'highlights','highlight'}))
  if isnumeric(cfg.highlight)
    plot_vector(x(cfg.highlight),y(cfg.highlight),'style',cfg.hlmarker,'Color',cfg.hlcolor,'markersize',cfg.hlmarkersize, 'linewidth',cfg.hllinewidth, 'markerfacecolor', cfg.hlfacecolor);
  else
    error('Unknown highlight type');
  end;
elseif strcmp(cfg.electrodes,'labels') || strcmp(cfg.showlabels,'yes')
  for i = 1:numChan
    plot_text(x(i), y(i), chanLabels{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', cfg.ecolor, 'FontSize', cfg.efontsize);
  end
elseif strcmp(cfg.electrodes,'numbers') || strcmp(cfg.showlabels,'numbers')
  for i = 1:numChan
    plot_text(x(i), y(i), int2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', cfg.ecolor, 'FontSize',cfg.efontsize);
  end
elseif strcmp(cfg.electrodes,'dotnum')
  for i = 1:numChan
    plot_text(x(i), y(i), int2str(i), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', cfg.ecolor, 'FontSize', cfg.efontsize);
  end
  if ischar(cfg.highlight)
    plot_vector(x, y, 'style', cfg.emarker, 'Color', cfg.ecolor, 'markersize', cfg.emarkersize);
  elseif isnumeric(cfg.highlight)
    normal = setdiff(1:length(x), cfg.highlight);
    plot_vector(x(normal), y(normal), 'style', cfg.emarker, 'Color', cfg.ecolor, 'markersize', cfg.emarkersize);
    plot_vector(x(cfg.highlight), y(cfg.highlight), 'style', cfg.hlmarker, 'Color', cfg.hlcolor, 'markersize', cfg.hlmarkersize, 'linewidth', cfg.hllinewidth, 'markerfacecolor', cfg.hlfacecolor);
  else
    error('Unknown highlight type');
  end;
end

if isfield(cfg.layout, 'outline')
  % plot the outline of the head, ears and nose
  for i=1:length(cfg.layout.outline)
    plot_vector(cfg.layout.outline{i}(:,1), cfg.layout.outline{i}(:,2), 'Color', cfg.headcolor, 'LineWidth', cfg.hlinewidth)
  end
end

% Write comment:
if isfield(cfg, 'comment')
  if strcmp(cfg.commentpos, 'title')
    title(cfg.comment, 'Fontsize', cfg.fontsize);
  else
    plot_text(x_COMNT, y_COMNT, cfg.comment, 'Fontsize', cfg.fontsize, 'HorizontalAlignment', HorAlign, 'VerticalAlignment', VerAlign);
  end
end

% plot colorbar:
if isfield(cfg, 'colorbar') && ~all(data == data(1))
  if strcmp(cfg.colorbar, 'yes')
    colorbar;
  elseif ~strcmp(cfg.colorbar, 'no')
    colorbar('location',cfg.colorbar);
  end
end

hold off
axis off
axis tight
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this allows to update the existing figure, this is to speed up realtime
% plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(cfg, 'update') && strcmp(cfg.update, 'yes')
  update.Xi = Xi;
  update.Yi = Yi;
  update.Zi = Zi;
  update.xi = xi;
  update.yi = yi;
  update.x  = x;
  update.y  = y;
  guidata(gcf, update);
end
