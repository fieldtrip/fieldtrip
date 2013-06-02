function [layout, cfg] = ft_prepare_layout(cfg, data)

% FT_PREPARE_LAYOUT loads or creates a 2-D layout of the channel locations. This layout
% is required for plotting the topographical distribution of the potential
% or field distribution, or for plotting timecourses in a topographical
% arrangement.
%
% Use as
%   layout = ft_prepare_layout(cfg, data)
%
% There are several ways in which a 2-D layout can be made: it can be read
% directly from a *.mat file containing a variable 'lay', it can be created
% based on 3-D electrode or gradiometer positions in the configuration or
% in the data, or it can be created based on the specification of an electrode
% or gradiometer file. Layouts can also come from an ASCII *.lay file, but
% this type of layout is no longer recommended.
%
% You can specify any one of the following configuration options
%   cfg.layout      filename containg the layout (.mat or .lay file)
%                   can also be a layout structure, which is simply returned
%                   as-is (see below for details)
%   cfg.rotate      number, rotation around the z-axis in degrees (default = [], which means automatic)
%   cfg.projection  string, 2D projection method can be 'stereographic', 'orthographic', 'polar', 'gnomic' or 'inverse' (default = 'polar')
%   cfg.elec        structure with electrode positions, or
%   cfg.elecfile    filename containing electrode positions
%   cfg.grad        structure with gradiometer definition, or
%   cfg.gradfile    filename containing gradiometer definition
%   cfg.output      filename to which the layout will be written (default = [])
%   cfg.montage     'no' or a montage structure (default = 'no')
%   cfg.image       filename, use an image to construct a layout (e.g. usefull for ECoG grids)
%   cfg.bw          if an image is used and bw = 1 transforms the image in black and white (default = 0, do not transform)
%   cfg.overlap     string, how to deal with overlapping channels when
%                   layout is constructed from a sensor configuration
%                   structure (can be 'shift' (shift the positions in 2D
%                   space to remove the overlap (default)), 'keep' (don't shift,
%                   retain the overlap), 'no' (throw error when overlap is
%                   present))
%  cfg.skipscale    'yes' or 'no', whether the scale should be included in the layout or not (default = 'no')
%  cfg.skipcomnt    'yes' or 'no', whether the comment should be included in the layout or not (default = 'no')
%
% Alternatively the layout can be constructed from either
%   data.elec     structure with electrode positions
%   data.grad     structure with gradiometer definition
%
% Alternatively you can specify the following layouts which will be
% generated for all channels present in the data. Note that these layouts
% are suitable for multiplotting, but not for topoplotting.
%   cfg.layout = 'ordered'  will give you a NxN ordered layout
%   cfg.layout = 'vertical' will give you a Nx1 ordered layout
%   cfg.layout = 'butterfly'  will give you a layout with all channels on top of each other
%
% The output layout structure will contain the following fields
%   layout.label   = Nx1 cell-array with channel labels
%   layout.pos     = Nx2 matrix with channel positions
%   layout.width   = Nx1 vector with the width of each box for multiplotting
%   layout.height  = Nx1 matrix with the height of each box for multiplotting
%   layout.mask    = optional cell-array with line segments that determine the area for topographic interpolation
%   layout.outline = optional cell-array with line segments that represent the head, nose, ears, sulci or other anatomical features
%
% See also FT_LAYOUTPLOT, FT_TOPOPLOTER, FT_TOPOPLOTTFR, FT_MULTIPLOTER, FT_MULTIPLOTTFR

% TODO switch to using planarchannelset function

% undocumented and non-recommended option (for SPM only)
%   cfg.style       string, '2d' or '3d' (default = '2d')

% Copyright (C) 2007-2009, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic check/initialization of input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
  data = [];
end
% ft_checkdata used to be called here in case data nargin>1, I moved this
% down to the branches of the big if-else-tree where data was actually
% used. speedup ~500ms (ES, dec2012)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default configuration options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg, 'rotate'),     cfg.rotate = [];                end  % [] => rotation is determined based on the type of sensors
if ~isfield(cfg, 'style'),      cfg.style = '2d';               end
if ~isfield(cfg, 'projection'), cfg.projection = 'polar';       end
if ~isfield(cfg, 'layout'),     cfg.layout = [];                end
if ~isfield(cfg, 'grad'),       cfg.grad = [];                  end
if ~isfield(cfg, 'elec'),       cfg.elec = [];                  end
if ~isfield(cfg, 'gradfile'),   cfg.gradfile = [];              end
if ~isfield(cfg, 'elecfile'),   cfg.elecfile = [];              end
if ~isfield(cfg, 'output'),     cfg.output = [];                end
if ~isfield(cfg, 'feedback'),   cfg.feedback = 'no';            end
if ~isfield(cfg, 'montage'),    cfg.montage = 'no';             end
if ~isfield(cfg, 'image'),      cfg.image = [];                 end
if ~isfield(cfg, 'bw'),         cfg.bw = 0;                     end
if ~isfield(cfg, 'channel'),    cfg.channel = 'all';            end
if ~isfield(cfg, 'skipscale'),  cfg.skipscale = 'no';           end
if ~isfield(cfg, 'skipcomnt'),  cfg.skipcomnt = 'no';           end
if ~isfield(cfg, 'overlap'),    cfg.overlap = 'shift';          end

cfg = ft_checkconfig(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to generate the layout structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

skipscale = strcmp(cfg.skipscale, 'yes'); % in general a scale is desired
skipcomnt = strcmp(cfg.skipcomnt, 'yes'); % in general a comment desired

if isa(cfg.layout, 'config')
  % convert the nested config-object back into a normal structure
  cfg.layout = struct(cfg.layout);
end

% ensure that there is a label field in the data, which is needed for
% ordered/vertical/butterfly modes
if nargin>1 && ~isfield(data, 'label') && isfield(data, 'labelcmb')
  data.label = unique(data.labelcmb(:));
end

% check whether cfg.layout already contains a valid layout structure (this can
% happen when higher level plotting functions are called with cfg.layout set to
% a layout structure)
if isstruct(cfg.layout) && isfield(cfg.layout, 'pos') && isfield(cfg.layout, 'label') && isfield(cfg.layout, 'width') && isfield(cfg.layout, 'height')
  layout = cfg.layout;
  
elseif isstruct(cfg.layout) && isfield(cfg.layout, 'pos') && isfield(cfg.layout, 'label') && (~isfield(cfg.layout, 'width') || ~isfield(cfg.layout, 'height'))
  layout = cfg.layout;
  % add width and height for multiplotting
  d = dist(layout.pos');
  nchans = length(layout.label);
  for i=1:nchans
    d(i,i) = inf; % exclude the diagonal
  end
  mindist = min(d(:));
  layout.width  = ones(nchans,1) * mindist * 0.8;
  layout.height = ones(nchans,1) * mindist * 0.6;
  
elseif isequal(cfg.layout, 'butterfly')
  if nargin>1 && ~isempty(data)
    % look at the data to determine the overlapping channels
    cfg.channel = ft_channelselection(cfg.channel, data.label);
    chanindx    = match_str(data.label, cfg.channel);
    nchan       = length(data.label(chanindx));
    layout.label   = data.label(chanindx);
  else
    nchan     = length(cfg.channel);
    layout.label = cfg.channel;
  end
  layout.pos     = zeros(nchan,2);  % centered at (0,0)
  layout.width   = ones(nchan,1) * 1.0;
  layout.height  = ones(nchan,1) * 1.0;
  layout.mask    = {};
  layout.outline = {};
  skipscale = true; % a scale is not desired
  skipcomnt = true; % a comment is initially not desired, or at least requires more thought
  
elseif isequal(cfg.layout, 'vertical')
  if nargin>1 && ~isempty(data)
    % look at the data to determine the overlapping channels
    data = ft_checkdata(data);
    cfg.channel = ft_channelselection(data.label, cfg.channel); % with this order order of channels stays the same
    [dum chanindx] = match_str(cfg.channel, data.label); % order of channels according to cfg specified by user
    nchan       = length(data.label(chanindx));
    layout.label   = data.label(chanindx);
  else
    nchan     = length(cfg.channel);
    layout.label = cfg.channel;
  end
  for i=1:(nchan+2)
    x = 0.5;
    y = 1-i/(nchan+1+2);
    layout.pos   (i,:) = [x y];
    layout.width (i,1) = 0.9;
    layout.height(i,1) = 0.9 * 1/(nchan+1+2);
    if i==(nchan+1)
      layout.label{i}   = 'SCALE';
    elseif i==(nchan+2)
      layout.label{i}   = 'COMNT';
    end
  end
  layout.mask    = {};
  layout.outline = {};
  
elseif isequal(cfg.layout, 'ordered')
  if nargin>1 && ~isempty(data)
    % look at the data to determine the overlapping channels
    data = ft_checkdata(data);
    cfg.channel = ft_channelselection(cfg.channel, data.label);
    chanindx    = match_str(data.label, cfg.channel);
    nchan       = length(data.label(chanindx));
    layout.label   = data.label(chanindx);
  else
    nchan     = length(cfg.channel);
    layout.label = cfg.channel;
  end
  ncol = ceil(sqrt(nchan))+1;
  nrow = ceil(sqrt(nchan))+1;
  k = 0;
  for i=1:nrow
    for j=1:ncol
      k = k+1;
      if k<=nchan
        x = (j-1)/ncol;
        y = (nrow-i-1)/nrow;
        layout.pos(k,:) = [x y];
        layout.width(k,1)  = 0.8 * 1/ncol;
        layout.height(k,1) = 0.8 * 1/nrow;
      end
    end
  end
  
  layout.label{end+1}  = 'SCALE';
  layout.width(end+1)  = 0.8 * 1/ncol;
  layout.height(end+1) = 0.8 * 1/nrow;
  x = (ncol-2)/ncol;
  y = 0/nrow;
  layout.pos(end+1,:) = [x y];
  
  layout.label{end+1}  = 'COMNT';
  layout.width(end+1)  = 0.8 * 1/ncol;
  layout.height(end+1) = 0.8 * 1/nrow;
  x = (ncol-1)/ncol;
  y = 0/nrow;
  layout.pos(end+1,:) = [x y];
  
  % try to generate layout from other configuration options
elseif ischar(cfg.layout)
  
  % layout file name specified
  
  if isempty(strfind(cfg.layout, '.'))
    
    cfg.layout = [cfg.layout '.mat'];
    if exist(cfg.layout, 'file')
      fprintf('layout file without .mat (or .lay) extension specified, appending .mat\n');
      layout = ft_prepare_layout(cfg);
    else
      cfg.layout = [cfg.layout(1:end-3) 'lay'];
      layout = ft_prepare_layout(cfg);
    end
    
  elseif ft_filetype(cfg.layout, 'matlab')
    
    fprintf('reading layout from file %s\n', cfg.layout);
    if ~exist(cfg.layout, 'file')
      error('the specified layout file %s was not found', cfg.layout);
    end
    tmp = load(cfg.layout, 'lay');
    if isfield(tmp, 'layout')
      layout = tmp.layout;
    elseif isfield(tmp, 'lay')
      layout = tmp.lay;
    else
      error('mat file does not contain a layout');
    end
    
  elseif ft_filetype(cfg.layout, 'layout')
    
    if exist(cfg.layout, 'file')
      fprintf('reading layout from file %s\n', cfg.layout);
      layout = readlay(cfg.layout);
    else
      warning_once(sprintf('layout file %s was not found on your path, attempting to use a similarly named .mat file instead',cfg.layout));
      cfg.layout = [cfg.layout(1:end-3) 'mat'];
      layout = ft_prepare_layout(cfg);
    end
    
  elseif ~ft_filetype(cfg.layout, 'layout')
    % assume that cfg.layout is an electrode file
    fprintf('creating layout from electrode file %s\n', cfg.layout);
    layout = sens2lay(ft_read_sens(cfg.layout), cfg.rotate, cfg.projection, cfg.style, cfg.overlap);
  end
  
elseif ischar(cfg.elecfile)
  fprintf('creating layout from electrode file %s\n', cfg.elecfile);
  layout = sens2lay(ft_read_sens(cfg.elecfile), cfg.rotate, cfg.projection, cfg.style, cfg.overlap);
  
elseif ~isempty(cfg.elec) && isstruct(cfg.elec)
  fprintf('creating layout from cfg.elec\n');
  layout = sens2lay(cfg.elec, cfg.rotate, cfg.projection, cfg.style, cfg.overlap);
  
elseif isfield(data, 'elec') && isstruct(data.elec)
  fprintf('creating layout from data.elec\n');
  data = ft_checkdata(data);
  layout = sens2lay(data.elec, cfg.rotate, cfg.projection, cfg.style, cfg.overlap);
  
elseif ischar(cfg.gradfile)
  fprintf('creating layout from gradiometer file %s\n', cfg.gradfile);
  layout = sens2lay(ft_read_sens(cfg.gradfile), cfg.rotate, cfg.projection, cfg.style, cfg.overlap);
  
elseif ~isempty(cfg.grad) && isstruct(cfg.grad)
  fprintf('creating layout from cfg.grad\n');
  layout = sens2lay(cfg.grad, cfg.rotate, cfg.projection, cfg.style, cfg.overlap);
  
elseif isfield(data, 'grad') && isstruct(data.grad)
  fprintf('creating layout from data.grad\n');
  data = ft_checkdata(data);
  layout = sens2lay(data.grad, cfg.rotate, cfg.projection, cfg.style, cfg.overlap);
  
elseif ~isempty(cfg.image) && isempty(cfg.layout)
  fprintf('reading background image from %s\n', cfg.image);
  [p,f,e] = fileparts(cfg.image);
  switch e
    case '.mat'
      tmp = whos('-file', cfg.image);
      load(cfg.image, tmp.name);
      eval(['img = ',tmp.name,';']);
    otherwise
      img = imread(cfg.image);
  end
  img = flipdim(img, 1); % in combination with "axis xy"
  
  figure
  bw = cfg.bw;
  
  if bw
    % convert to greyscale image
    img = mean(img, 3);
    imagesc(img);
    colormap gray
  else
    % plot as RGB image
    image(img);
  end
  
  hold on
  axis equal
  axis off
  axis xy
  
  % get the electrode positions
  pos = zeros(0,2);
  electrodehelp = [ ...
    '-----------------------------------------------------\n' ...
    'specify electrode locations\n' ...
    'press the right mouse button to add another electrode\n' ...
    'press backspace on the keyboard to remove the last electrode\n' ...
    'press "q" on the keyboard to continue\n' ...
    ];
  again = 1;
  while again
    fprintf(electrodehelp)
    disp(round(pos)); % values are integers/pixels
    try
      [x, y, k] = ginput(1);
    catch
      % this happens if the figure is closed
      return;
    end
    
    switch k
      case 1
        pos = cat(1, pos, [x y]);
        % add it to the figure
        plot(x, y, 'b.');
        plot(x, y, 'yo');
        
      case 8
        if size(pos,1)>0
          % remove the last point
          pos = pos(1:end-1,:);
          % completely redraw the figure
          cla
          h = image(img);
          hold on
          axis equal
          axis off
          plot(pos(:,1), pos(:,2), 'b.');
          plot(pos(:,1), pos(:,2), 'yo');
        end
        
      case 'q'
        again = 0;
        
      otherwise
        warning('invalid button (%d)', k);
    end
  end
  
  % get the mask outline
  polygon = {};
  thispolygon = 1;
  polygon{thispolygon} = zeros(0,2);
  maskhelp = [ ...
    '------------------------------------------------------------------------\n' ...
    'specify polygons for masking the topgraphic interpolation\n' ...
    'press the right mouse button to add another point to the current polygon\n' ...
    'press backspace on the keyboard to remove the last point\n' ...
    'press "c" on the keyboard to close this polygon and start with another\n' ...
    'press "q" on the keyboard to continue\n' ...
    ];
  again = 1;
  while again
    fprintf(maskhelp);
    fprintf('\n');
    for i=1:length(polygon)
      fprintf('polygon %d has %d points\n', i, size(polygon{i},1));
    end
    
    try
      [x, y, k] = ginput(1);
    catch
      % this happens if the figure is closed
      return;
    end
    
    switch k
      case 1
        polygon{thispolygon} = cat(1, polygon{thispolygon}, [x y]);
        % add the last line segment to the figure
        if size(polygon{thispolygon},1)>1
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
        end
        plot(x, y, 'g.-');
        
      case 8 % backspace
        if size(polygon{thispolygon},1)>0
          % remove the last point
          polygon{thispolygon} = polygon{thispolygon}(1:end-1,:);
          % completely redraw the figure
          cla
          h = image(img);
          hold on
          axis equal
          axis off
          % plot the electrode positions
          plot(pos(:,1), pos(:,2), 'b.');
          plot(pos(:,1), pos(:,2), 'yo');
          for i=1:length(polygon)
            x = polygon{i}(:,1);
            y = polygon{i}(:,2);
            if i~=thispolygon
              % close the polygon in the figure
              x(end) = x(1);
              y(end) = y(1);
            end
            plot(x, y, 'g.-');
          end
        end
        
      case 'c'
        if size(polygon{thispolygon},1)>0
          % close the polygon
          polygon{thispolygon}(end+1,:) = polygon{thispolygon}(1,:);
          % close the polygon in the figure
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
          plot(x, y, 'g.-');
          % switch to the next polygon
          thispolygon = thispolygon + 1;
          polygon{thispolygon} = zeros(0,2);
        end
        
      case 'q'
        if size(polygon{thispolygon},1)>0
          % close the polygon
          polygon{thispolygon}(end+1,:) = polygon{thispolygon}(1,:);
          % close the polygon in the figure
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
          plot(x, y, 'g.-');
        end
        again = 0;
        
      otherwise
        warning('invalid button (%d)', k);
    end
  end
  % remember this set of polygons as the mask
  mask = polygon;
  
  
  % get the outline, e.g. head shape and sulci
  polygon = {};
  thispolygon = 1;
  polygon{thispolygon} = zeros(0,2);
  maskhelp = [ ...
    '-----------------------------------------------------------------------------------\n' ...
    'specify polygons for adding outlines (e.g. head shape and sulci) to the layout\n' ...
    'press the right mouse button to add another point to the current polygon\n' ...
    'press backspace on the keyboard to remove the last point\n' ...
    'press "c" on the keyboard to close this polygon and start with another\n' ...
    'press "n" on the keyboard to start with another without closing the current polygon\n' ...
    'press "q" on the keyboard to continue\n' ...
    ];
  again = 1;
  while again
    fprintf(maskhelp);
    fprintf('\n');
    for i=1:length(polygon)
      fprintf('polygon %d has %d points\n', i, size(polygon{i},1));
    end
    
    try
      [x, y, k] = ginput(1);
    catch
      % this happens if the figure is closed
      return;
    end
    
    switch k
      case 1
        polygon{thispolygon} = cat(1, polygon{thispolygon}, [x y]);
        % add the last line segment to the figure
        if size(polygon{thispolygon},1)>1
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
        end
        plot(x, y, 'm.-');
        
      case 8 % backspace
        if size(polygon{thispolygon},1)>0
          % remove the last point
          polygon{thispolygon} = polygon{thispolygon}(1:end-1,:);
          % completely redraw the figure
          cla
          h = image(img);
          hold on
          axis equal
          axis off
          % plot the electrode positions
          plot(pos(:,1), pos(:,2), 'b.');
          plot(pos(:,1), pos(:,2), 'yo');
          for i=1:length(polygon)
            x = polygon{i}(:,1);
            y = polygon{i}(:,2);
            if i~=thispolygon
              % close the polygon in the figure
              x(end) = x(1);
              y(end) = y(1);
            end
            plot(x, y, 'm.-');
          end
        end
        
      case 'c'
        if size(polygon{thispolygon},1)>0
          x = polygon{thispolygon}(1,1);
          y = polygon{thispolygon}(1,2);
          polygon{thispolygon} = cat(1, polygon{thispolygon}, [x y]);
          % add the last line segment to the figure
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
          plot(x, y, 'm.-');
          % switch to the next polygon
          thispolygon = thispolygon + 1;
          polygon{thispolygon} = zeros(0,2);
        end
        
      case 'n'
        if size(polygon{thispolygon},1)>0
          % switch to the next polygon
          thispolygon = thispolygon + 1;
          polygon{thispolygon} = zeros(0,2);
        end
        
      case 'q'
        again = 0;
        
      otherwise
        warning('invalid button (%d)', k);
    end
  end
  % remember this set of polygons as the outline
  outline = polygon;
  
  % convert electrode positions into a layout structure
  layout.pos = pos;
  nchans = size(pos,1);
  for i=1:nchans
    layout.label{i,1} = sprintf('chan%03d', i);
  end
  % add width and height for multiplotting
  d = dist(pos');
  for i=1:nchans
    d(i,i) = inf; % exclude the diagonal
  end
  mindist = min(d(:));
  layout.width  = ones(nchans,1) * mindist * 0.8;
  layout.height = ones(nchans,1) * mindist * 0.6;
  % add mask and outline polygons
  layout.mask    = mask;
  layout.outline = outline;
  
else
  error('no layout detected, please specify cfg.layout')
end

% FIXME there is a conflict between the use of cfg.style here and in topoplot
if ~strcmp(cfg.style, '3d')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % check whether outline and mask are available
  % if not, add default "circle with triangle" to resemble the head
  % in case of "circle with triangle", the electrode positions should also be
  % scaled
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isfield(layout, 'outline') || ~isfield(layout, 'mask')
    rmax  = 0.5;
    l     = 0:2*pi/100:2*pi;
    HeadX = cos(l).*rmax;
    HeadY = sin(l).*rmax;
    NoseX = [0.18*rmax 0 -0.18*rmax];
    NoseY = [rmax-.004 rmax*1.15 rmax-.004];
    EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
    EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
    % Scale the electrode positions to fit within a unit circle, i.e. electrode radius = 0.45
    ind_scale = strmatch('SCALE', layout.label);
    ind_comnt = strmatch('COMNT', layout.label);
    sel = setdiff(1:length(layout.label), [ind_scale ind_comnt]); % these are excluded for scaling
    x = layout.pos(sel,1);
    y = layout.pos(sel,2);
    xrange = range(x);
    yrange = range(y);
    % First scale the width and height of the box for multiplotting
    layout.width  = layout.width./xrange;
    layout.height = layout.height./yrange;
    % Then shift and scale the electrode positions
    layout.pos(:,1) = 0.9*((layout.pos(:,1)-min(x))/xrange-0.5);
    layout.pos(:,2) = 0.9*((layout.pos(:,2)-min(y))/yrange-0.5);
    % Define the outline of the head, ears and nose
    layout.outline{1} = [HeadX(:) HeadY(:)];
    layout.outline{2} = [NoseX(:) NoseY(:)];
    layout.outline{3} = [ EarX(:)  EarY(:)];
    layout.outline{4} = [-EarX(:)  EarY(:)];
    % Define the anatomical mask based on a circular head
    layout.mask{1} = [HeadX(:) HeadY(:)];
  end
end


% make the subset as specified in cfg.channel
cfg.channel = ft_channelselection(cfg.channel, setdiff(layout.label, {'COMNT', 'SCALE'}));  % COMNT and SCALE are not really channels
chansel = match_str(layout.label, cat(1, cfg.channel(:), 'COMNT', 'SCALE'));                % include COMNT and SCALE, keep all channels in the order of the layout
% return the layout for the subset of channels
layout.pos    = layout.pos(chansel,:);
layout.label  = layout.label(chansel);
if ~strcmp(cfg.style, '3d')
  % these don't exist for the 3D layout
  layout.width  = layout.width(chansel);
  layout.height = layout.height(chansel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the montage, i.e. combine bipolar channels into a new representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(cfg.montage, 'no')
  Norg = length(cfg.montage.labelorg);
  Nnew = length(cfg.montage.labelnew);
  
  for i=1:Nnew
    cfg.montage.tra(i,:) = abs(cfg.montage.tra(i,:));
    cfg.montage.tra(i,:) = cfg.montage.tra(i,:) ./ sum(cfg.montage.tra(i,:));
  end
  % pretend it is a sensor structure, this achieves averaging after channel matching
  tmp.tra   = layout.pos;
  tmp.label = layout.label;
  new = ft_apply_montage(tmp, cfg.montage);
  layout.pos   = new.tra;
  layout.label = new.label;
  % do the same for the width and height
  tmp.tra = layout.width(:);
  new = ft_apply_montage(tmp, cfg.montage);
  layout.width = new.tra;
  tmp.tra = layout.height(:);
  new = ft_apply_montage(tmp, cfg.montage);
  layout.height = new.tra;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add axes positions for comments and scale information if required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~any(strcmp('COMNT', layout.label)) && strcmpi(cfg.style, '2d') && ~skipcomnt
  % add a placeholder for the comment in the upper left corner
  layout.label{end+1}  = 'COMNT';
  layout.width(end+1)  = mean(layout.width);
  layout.height(end+1) = mean(layout.height);
  X                 = min(layout.pos(:,1));
  Y                 = max(layout.pos(:,2));
  Y                 = min(layout.pos(:,2));
  layout.pos(end+1,:)  = [X Y];
elseif any(strcmp('COMNT', layout.label)) && skipcomnt
  % remove the scale entry
  sel = find(strcmp('COMNT', layout.label));
  layout.label(sel) = [];
  layout.pos(sel,:) = [];
  layout.width(sel) = [];
  layout.height(sel) = [];
end

if ~any(strcmp('SCALE', layout.label)) && strcmpi(cfg.style, '2d') && ~skipscale
  % add a placeholder for the scale in the upper right corner
  layout.label{end+1}  = 'SCALE';
  layout.width(end+1)  = mean(layout.width);
  layout.height(end+1) = mean(layout.height);
  X                 = max(layout.pos(:,1));
  Y                 = max(layout.pos(:,2));
  Y                 = min(layout.pos(:,2));
  layout.pos(end+1,:)  = [X Y];
elseif any(strcmp('SCALE', layout.label)) && skipscale
  % remove the scale entry
  sel = find(strcmp('SCALE', layout.label));
  layout.label(sel) = [];
  layout.pos(sel,:) = [];
  layout.width(sel) = [];
  layout.height(sel) = [];
end

% to plot the layout for debugging, you can use this code snippet
if strcmp(cfg.feedback, 'yes') && strcmpi(cfg.style, '2d')
  tmpcfg = [];
  tmpcfg.layout = layout;
  ft_layoutplot(tmpcfg);
end

% to write the layout to a text file, you can use this code snippet
if ~isempty(cfg.output) && strcmpi(cfg.style, '2d')
  fprintf('writing layout to ''%s''\n', cfg.output);
  fid = fopen(cfg.output, 'wt');
  for i=1:numel(layout.label)
    fprintf(fid, '%d %f %f %f %f %s\n', i, layout.pos(i,1), layout.pos(i,2), layout.width(i), layout.height(i), layout.label{i});
  end
  fclose(fid);
elseif ~isempty(cfg.output) && strcmpi(cfg.style, '3d')
  % the layout file format does not support 3D positions, furthermore for
  % a 3D layout the width and height are currently set to NaN
  error('writing a 3D layout to an output file is not supported');
end

% ensure proper format of some of label (see bug 1909 -roevdmei)
layout.label  = layout.label(:);


% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble history layout

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% read the layout information from the ascii file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = readlay(filename)
if ~exist(filename, 'file')
  error(sprintf('could not open layout file: %s', filename));
end
[chNum,X,Y,Width,Height,Lbl,Rem] = textread(filename,'%f %f %f %f %f %q %q');

if length(Rem)<length(Lbl)
  Rem{length(Lbl)} = [];
end

for i=1:length(Lbl)
  if ~isempty(Rem{i})
    % this ensures that channel names with a space in them are also supported (i.e. Neuromag)
    Lbl{i} = [Lbl{i} ' ' Rem{i}];
  end
end
layout.pos    = [X Y];
layout.width  = Width;
layout.height = Height;
layout.label  = Lbl;
return % function readlay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% convert 3D electrode positions into 2D layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = sens2lay(sens, rz, method, style, overlap)

% remove the balancing from the sensor definition, e.g. 3rd order gradients, PCA-cleaned data or ICA projections
% this should not be necessary anymore, because the sensor description is
% up-to-date, i.e. explicit information with respect to the channel
% positions is present
%sens = undobalancing(sens);

fprintf('creating layout for %s system\n', ft_senstype(sens));

% apply rotation
if isempty(rz)
  switch ft_senstype(sens)
    case {'ctf151', 'ctf275', 'bti148', 'bti248', 'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar', 'yokogawa160', 'yokogawa160_planar', 'yokogawa64', 'yokogawa64_planar', 'yokogawa440', 'yokogawa440_planar', 'magnetometer', 'meg'}
      rz = 90;
    case {'neuromag122', 'neuromag306'}
      rz = 0;
    case 'electrode'
      rz = 90;
    otherwise
      rz = 0;
  end
end
sens.chanpos = warp_apply(rotate([0 0 rz]), sens.chanpos, 'homogenous');

% determine the 3D channel positions
pnt   = sens.chanpos;
label = sens.label;

if strcmpi(style, '3d')
  layout.pos   = pnt;
  layout.label = label;
else
  prj = elproj(pnt, method);
  
  % this copy will be used to determine the minimum distance between channels
  % we need a copy because prj retains the original positions, and
  % prjForDist might need to be changed if the user wants to keep
  % overlapping channels
  prjForDist = prj;
  
  % check whether many channels occupy identical positions, if so shift
  % them around if requested
  if size(unique(prj,'rows'),1) / size(prj,1) < 0.8
    if strcmp(overlap, 'shift')
      warning_once('the specified sensor configuration has many overlapping channels, creating a layout by shifting them around (use a template layout for better control over the positioning)');
      prj = shiftxy(prj', 0.2)';
      prjForDist = prj;
    elseif strcmp(overlap, 'no')
      error('the specified sensor configuration has many overlapping channels, you specified not to allow that');
    elseif strcmp(overlap, 'keep')
      prjForDist = unique(prj, 'rows');
    else
      error('unknown value for cfg.overlap = ''%s''', overlap);
    end
  end
  
  d = dist(prjForDist');
  d(logical(eye(size(d)))) = inf;
  
  % This is a fix for .sfp files, containing positions of 'fiducial
  % electrodes'. Their presence determines the minimum distance between
  % projected electrode positions, leading to very small boxes.
  % This problem has been detected and reported by Matt Mollison.
  % FIXME: consider changing the box-size being determined by mindist
  % by a mean distance or so; this leads to overlapping boxes, but that
  % also exists for some .lay files
  if any(strmatch('Fid', label))
    tmpsel = strmatch('Fid', label);
    d(tmpsel, :) = inf;
    d(:, tmpsel) = inf;
  end
  
  % take mindist as the median of the first quartile of closest channel pairs with non-zero distance
  mindist = min(d); % get closest neighbour for all channels
  mindist = sort(mindist(mindist>1e-6),'ascend');
  mindist = mindist(1:round(numel(label)/4));
  mindist = median(mindist);
  
  X = prj(:,1);
  Y = prj(:,2);
  Width  = ones(size(X)) * mindist * 0.8;
  Height = ones(size(X)) * mindist * 0.6;
  layout.pos    = [X Y];
  layout.width  = Width;
  layout.height = Height;
  layout.label  = label;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% shift 2D positions around so that the minimum distance between any pair
% is mindist
%
% Credit for this code goes to Laurence Hunt at UCL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xy = shiftxy(xy,mindist)

x = xy(1,:);
y = xy(2,:);

l=1;
i=1; %filler
mindist = mindist/0.999; % limits the number of loops
while (~isempty(i) && l<50)
  xdiff = repmat(x,length(x),1) - repmat(x',1,length(x));
  ydiff = repmat(y,length(y),1) - repmat(y',1,length(y));
  xydist= sqrt(xdiff.^2 + ydiff.^2); %euclidean distance between all sensor pairs
  
  [i,j] = find(xydist<mindist*0.999);
  rm=(i<=j); i(rm)=[]; j(rm)=[]; %only look at i>j
  
  for m = 1:length(i);
    if (xydist(i(m),j(m)) == 0)
      diffvec = [mindist./sqrt(2) mindist./sqrt(2)];
    else
      xydiff = [xdiff(i(m),j(m)) ydiff(i(m),j(m))];
      diffvec = xydiff.*mindist./xydist(i(m),j(m)) - xydiff;
    end
    x(i(m)) = x(i(m)) - diffvec(1)/2;
    y(i(m)) = y(i(m)) - diffvec(2)/2;
    x(j(m)) = x(j(m)) + diffvec(1)/2;
    y(j(m)) = y(j(m)) + diffvec(2)/2;
  end
  l = l+1;
end

xy = [x; y];

