function [layout, cfg] = ft_prepare_layout(cfg, data)

% FT_PREPARE_LAYOUT loads or creates a 2-D layout of the channel locations. This
% layout is required for plotting the topographical distribution of the potential or
% field distribution, or for plotting timecourses in a topographical arrangement.
%
% Use as
%   layout = ft_prepare_layout(cfg)
% or
%   layout = ft_prepare_layout(cfg, data)
% where the optional data input argument is any of the FieldTrip data structures.
%
% This returns a layout structure with the following elements
%   layout.pos     = Nx2 matrix with the position where each channel should be plotted
%   layout.label   = Nx1 cell-array with the channel labels
%   layout.width   = Nx1 vector with the width of each box for multiplotting
%   layout.height  = Nx1 vector with the height of each box for multiplotting
%   layout.mask    = optional cell-array with line segments that determine the area for topographic interpolation
%   layout.outline = optional cell-array with line segments that represent the head, nose, ears, sulci or other anatomical features
%   layout.color   = optional Nx3 matrix with rgb values for the channels' color, for fine-grained color behavior
%
% There are several ways in which a 2-D layout can be made:
% 1) it can be read directly from a layout file
% 2) it can be created on basis of an image or photo,
% 3) it can be created from a projection of the 3-D sensor positions in the data, in the configuration, or in an electrode, gradiometer or optode file.
%
% Layout files are MATLAB *.mat files containing a single structure representing the layout
% (see above). The layout file can also be an ASCII file with the extension *.lay, although
% this file format is no longer recommended, since there is less control over the outline
% of the head and the mask within which the interpolation is done. A large number of
% template layout files is provided in the fieldtrip/template/layout directory. See
% also http://www.fieldtriptoolbox.org/template/layout
%
% You can specify any one of the following configuration options
%   cfg.layout      = filename containg the input layout (*.mat or *.lay file), this can also be a layout
%                     structure, which is simply returned as-is (see below for details)
%   cfg.output      = filename (ending in .mat or .lay) to which the layout will be written (default = [])
%   cfg.feedback    = 'yes' or 'no', whether to show an image of the layout (default = 'no')
%   cfg.elec        = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad        = structure with gradiometer definition or filename, see FT_READ_SENS
%   cfg.opto        = structure with optode definition or filename, see FT_READ_SENS
%   cfg.rotate      = number, rotation around the z-axis in degrees (default = [], which means automatic)
%   cfg.center      = string, center and scale the electrodes in the sphere that represents the head, can be 'yes' or 'no' (default = 'no')
%   cfg.projection  = string, 2D projection method can be 'stereographic', 'orthographic', 'polar' or 'gnomic' (default = 'polar')
%                     When 'orthographic', cfg.viewpoint can be used to indicate to specificy projection (keep empty for legacy projection)
%   cfg.viewpoint   = string indicating the view point that is used for orthographic projection of 3-D sensor
%                     positions to the 2-D plane. The possible viewpoints are
%                     'left'      - left  sagittal view,     L=anterior, R=posterior, top=top, bottom=bottom
%                     'right'     - right sagittal view,     L=posterior, R=anterior, top=top, bottom=bottom
%                     'topleft'   - view from the top top,   L=anterior, R=posterior, top=top, bottom=bottom
%                     'topright'  - view from the top right, L=posterior, R=anterior, top=top, bottom=bottom
%                     'inferior'  - inferior axial view,     L=R, R=L, top=anterior, bottom=posterior
%                     'superior'  - superior axial view,     L=L, R=R, top=anterior, bottom=posterior
%                     'anterior'  - anterior  coronal view,  L=R, R=L, top=top, bottom=bottom
%                     'posterior' - posterior coronal view,  L=L, R=R, top=top, bottom=bottom
%                     'auto'      - automatic guess of the most optimal of the above
%                      tip: use cfg.viewpoint = 'auto' per iEEG electrode grid/strip/depth for more accurate results
%                      tip: to obtain an overview of all iEEG electrodes, choose superior/inferior, use cfg.headshape/mri, and plot using FT_LAYOUTPLOT with cfg.box/mask = 'no'
%   cfg.outline     = string, how to create the outline, can be 'circle', 'doublecirclecross', 'helmet', 'square', 'convex', 'headshape', 'mri' or 'no' (default is automatic)
%   cfg.mask        = string, how to create the mask, can be 'circle', 'extended', 'square', 'convex', 'headshape', 'mri' or 'no' (default is automatic)
%   cfg.headshape   = surface mesh (for example pial or head) to be used for generating an outline, see FT_READ_HEADSHAPE for details
%   cfg.mri         = segmented anatomical MRI to be used for generating an outline, see FT_READ_MRI and FT_VOLUMESEGMENT for details
%   cfg.montage     = 'no' or a montage structure (default = 'no')
%   cfg.image       = filename, use an image to construct a layout (useful for ECoG grids)
%   cfg.bw          = 'yes' or 'no', if an image is used and this option is true, the image is transformed in black and white (default = 'no', i.e. do not transform)
%   cfg.overlap     = string, how to deal with overlapping channels when the layout is constructed from a sensor configuration structure. This can be
%                     'shift'  - shift the positions in 2D space to remove the overlap (default)
%                     'keep'   - do not shift, retain the overlap
%                     'no'     - throw an error when overlap is present
%   cfg.channel     = 'all', or Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%   cfg.boxchannel  = 'all', or Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%                      specificies channels to use for determining channel box size (default = 'all', recommended for MEG/EEG, a selection is recommended for iEEG)
%   cfg.skipscale   = 'yes' or 'no', whether the scale should be included in the layout or not (default = 'no')
%   cfg.skipcomnt   = 'yes' or 'no', whether the comment should be included in the layout or not (default = 'no')
%   cfg.color       = empty, 'spatial', or Nx3 matrix, if non-empty, an Nx3 color matrix based on the position 
%                     of the sensors will be added (default = [])
%
% If you use cfg.headshape or cfg.mri to create a headshape outline, the input
% geometry should be expressed in the same units and coordinate system as the input
% sensors.
%
% Alternatively the layout can be constructed from either one of these in the input data structure:
%   data.elec     = structure with electrode positions
%   data.grad     = structure with gradiometer definition
%   data.opto     = structure with optode definition
%
% Alternatively you can specify the following options for systematic layouts which
% will be generated for all channels present in the data. Note that these layouts are
% only suitable for multiplotting, not for topoplotting.
%   cfg.layout = 'ordered'    will give you a NxN ordered layout
%   cfg.layout = 'vertical'   will give you a Nx1 ordered layout
%   cfg.layout = 'horizontal' will give you a 1xN ordered layout
%   cfg.layout = 'butterfly'  will give you a layout with all channels on top of each other
%   cfg.layout = 'circular'   will distribute the channels on a circle
%   cfg.width  = scalar (default is automatic)
%   cfg.height = scalar (default is automatic)
%
% For an sEEG shaft the option cfg.layout='vertical' or 'horizontal' is useful to
% represent the channels in a linear sequence . In this case you can also specify the
% direction of the shaft as going from left-to-right, top-to-bottom, etc.
%   cfg.direction = string, can be any of 'LR', 'RL' (for horizontal), 'TB', 'BT' (for vertical)
%
% For an ECoG grid the option cfg.layout='ordered' is useful to represent the
% channels in a grid array. In this case you can also specify the number of rows
% and/or columns and hwo the channels increment over the grid (e.g. first
% left-to-right, then top-to-bottom). You can check the channel order of your grid
% using FT_PLOT_LAYOUT.
%   cfg.rows      = number of rows (default is automatic)
%   cfg.columns   = number of columns (default is automatic)
%   cfg.direction = string, can be any of 'LRTB', 'RLTB', 'LRBT', 'RLBT', 'TBLR', 'TBRL', 'BTLR', 'BTRL' (default = 'LRTB')
%
% See also FT_TOPOPLOTER, FT_TOPOPLOTTFR, FT_MULTIPLOTER, FT_MULTIPLOTTFR, FT_PLOT_LAYOUT

% undocumented and non-recommended option (for SPM only)
%   cfg.style       string, '2d' or '3d' (default = '2d')

% Copyright (C) 2007-2020, Robert Oostenveld
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
ft_preamble loadvar    data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the data can be passed as input argument or can be read from disk
hasdata = exist('data', 'var') && ~isempty(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic check/initialization of input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~hasdata
  data = struct([]);
else
  % check if the input data is valid for this function
  data = ft_checkdata(data);
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamed',    {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed',    {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed',    {'optofile', 'opto'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default configuration options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.rotate       = ft_getopt(cfg, 'rotate',     []); % [] => default rotation is determined based on the type of sensors
cfg.center       = ft_getopt(cfg, 'center',     'no');
cfg.style        = ft_getopt(cfg, 'style',      '2d');
cfg.projection   = ft_getopt(cfg, 'projection', 'polar');
cfg.layout       = ft_getopt(cfg, 'layout',     []);
cfg.grad         = ft_getopt(cfg, 'grad',       []);
cfg.elec         = ft_getopt(cfg, 'elec',       []);
cfg.opto         = ft_getopt(cfg, 'opto',       []);
cfg.output       = ft_getopt(cfg, 'output',     []);
cfg.feedback     = ft_getopt(cfg, 'feedback',   'no');
cfg.montage      = ft_getopt(cfg, 'montage',    'no');
cfg.image        = ft_getopt(cfg, 'image',      []);
cfg.mesh         = ft_getopt(cfg, 'mesh',       []); % experimental, should only work with meshes defined in 2D
cfg.bw           = ft_getopt(cfg, 'bw',         'no');
cfg.channel      = ft_getopt(cfg, 'channel',    'all');
cfg.skipscale    = ft_getopt(cfg, 'skipscale',  []); % see below
cfg.skipcomnt    = ft_getopt(cfg, 'skipcomnt',  []); % see below
cfg.boxchannel   = ft_getopt(cfg, 'boxchannel', cfg.channel);
cfg.overlap      = ft_getopt(cfg, 'overlap',    'shift');
cfg.viewpoint    = ft_getopt(cfg, 'viewpoint',  []);
cfg.headshape    = ft_getopt(cfg, 'headshape',  []); % separate form cfg.mesh
cfg.mri          = ft_getopt(cfg, 'mri',        []);
cfg.outline      = ft_getopt(cfg, 'outline',    []); % default is handled below
cfg.mask         = ft_getopt(cfg, 'mask',       []); % default is handled below
cfg.width        = ft_getopt(cfg, 'width',      []);
cfg.height       = ft_getopt(cfg, 'height',     []);
cfg.commentpos   = ft_getopt(cfg, 'commentpos', 'layout');
cfg.scalepos     = ft_getopt(cfg, 'scalepos',   'layout');
cfg.color        = ft_getopt(cfg, 'color',      []);

if isempty(cfg.skipscale)
  if ischar(cfg.layout) && any(strcmp(cfg.layout, {'ordered', 'vertical', 'horizontal', 'butterfly', 'circular'}))
    cfg.skipscale = 'yes';
  else
    cfg.skipscale = 'no';
  end
end

if isempty(cfg.skipcomnt)
  if ischar(cfg.layout) && any(strcmp(cfg.layout, {'ordered', 'vertical', 'horizontal', 'butterfly', 'circular'}))
    cfg.skipcomnt = 'yes';
  else
    cfg.skipcomnt = 'no';
  end
end

if isempty(cfg.outline) || strcmp(cfg.outline, 'yes')
  % determine the most suitable shape of the outline
  if ~isempty(cfg.headshape)
    cfg.outline = 'headshape';
  elseif ~isempty(cfg.mri)
    cfg.outline = 'mri';
  elseif ischar(cfg.layout) && any(strcmp(cfg.layout, {'horizontal', 'vertical', 'ordered'}))
    cfg.outline = 'no';
  elseif ~strcmp(cfg.projection, 'orthographic')
    cfg.outline = 'circle';
  else
    cfg.outline = 'no';
  end
end

if isempty(cfg.mask) || strcmp(cfg.mask, 'yes')
  % determine the most suitable shape of the mask
  if ~isempty(cfg.headshape)
    cfg.mask = 'headshape';
  elseif ~isempty(cfg.mri)
    cfg.mask = 'mri';
  elseif ischar(cfg.layout) && any(strcmp(cfg.layout, {'horizontal', 'vertical', 'ordered'}))
    cfg.mask = 'square';
  elseif ~strcmp(cfg.projection, 'orthographic')
    cfg.mask = 'circle';
  else
    cfg.mask = 'convex';
  end
end

% headshape/mri are mutually exclusive
if ~isempty(cfg.headshape) && ~isempty(cfg.mri)
  ft_error('cfg.headshape and cfg.mri are mutually exclusive, please use only one of the two')
end
% cfg.viewpoint can only be used together with cfg.projection = 'orthographic'
if ~isempty(cfg.viewpoint) && ~isequal(cfg.projection, 'orthographic')
  ft_error('cfg.viewpoint can only be used in the case of orthographic projection')
end
if ~isempty(cfg.viewpoint) && ~isempty(cfg.rotate)
  ft_error('cfg.viewpoint and cfg.rotate are mutually exclusive, please use only one of the two')
end

% update the selection of channels according to the data
if hasdata && isfield(data, 'topolabel')
  cfg.channel = ft_channelselection(cfg.channel, data.topolabel);
elseif hasdata && isfield(data, 'label')
  cfg.channel = ft_channelselection(cfg.channel, data.label);
elseif hasdata && isfield(data, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
end

if ischar(cfg.layout) && strcmp(cfg.layout, 'vertical')
  cfg.direction = ft_getopt(cfg, 'direction', 'TB'); % default is top-to-bottom
elseif ischar(cfg.layout) && strcmp(cfg.layout, 'horizontal')
  cfg.direction = ft_getopt(cfg, 'direction', 'LR'); % default is left-to-right
elseif ischar(cfg.layout) && strcmp(cfg.layout, 'ordered')
  cfg.direction = ft_getopt(cfg, 'direction', 'LRTB'); % default is left-to-right, then top-to-bottom
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to generate the layout structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

skipscale = istrue(cfg.skipscale); % in general a scale is desired
skipcomnt = istrue(cfg.skipcomnt); % in general a comment desired

% ensure that there is a label field in the data, which is needed for
% ordered/vertical//horizontal/butterfly modes
if hasdata && ~isfield(data, 'label') && isfield(data, 'labelcmb')
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
  
elseif isequal(cfg.layout, 'circular')
  rho = ft_getopt(cfg, 'rho', []);
  
  if hasdata && ~isempty(data)
    % look at the data to determine the overlapping channels
    cfg.channel  = ft_channelselection(cfg.channel, data.label);
    chanindx     = match_str(data.label, cfg.channel);
    nchan        = length(data.label(chanindx));
    layout.label = data.label(chanindx);
  else
    assert(iscell(cfg.channel), 'cfg.channel should be a valid set of channels');
    nchan        = length(cfg.channel);
    layout.label = cfg.channel;
  end
  
  if isempty(rho)
    % do an equally spaced layout, starting at 12 o'clock, going clockwise
    rho = linspace(0,1,nchan+1);
    rho = 2.*pi.*rho(1:end-1);
  else
    if numel(rho) ~= nchan
      ft_error('the number of elements in the polar angle vector should be equal to the number of channels');
    end
    
    % convert to radians
    rho = 2.*pi.*rho./360;
  end
  x   = sin(rho);
  y   = cos(rho);
  
  layout.pos     = [x(:) y(:)];
  layout.width   = ones(nchan,1) * 0.01;
  layout.height  = ones(nchan,1) * 0.01;
  layout.mask    = {};
  layout.outline = {};
  
elseif isequal(cfg.layout, 'butterfly')
  if hasdata && ~isempty(data)
    % look at the data to determine the channels to be plotted
    cfg.channel  = ft_channelselection(cfg.channel, data.label);
    chanindx     = match_str(data.label, cfg.channel);
    nchan        = length(data.label(chanindx));
    layout.label = data.label(chanindx);
  else
    assert(iscell(cfg.channel), 'cfg.channel should be a valid set of channels');
    nchan        = length(cfg.channel);
    layout.label = cfg.channel;
  end
  layout.pos     = zeros(nchan,2);  % centered at (0,0)
  layout.width   = ones(nchan,1) * 1.0;
  layout.height  = ones(nchan,1) * 1.0;
  layout.mask    = {};
  layout.outline = {};
  
  if ~istrue(cfg.skipscale)
    layout.label{end+1}  = 'SCALE';
    layout.pos(end+1, :) = layout.pos(1,:);
    layout.width(end+1)  = layout.width(1);
    layout.height(end+1) = layout.height(1);
  end
  
  if strcmp(cfg.color, 'spatial')
    try
      % make one with the default settings to copy the spatial colors
      tmpcfg = keepfields(cfg, {'channel', 'color', 'skipcomnt', 'skipscale'});
      tmpcfg = ft_setopt(tmpcfg, 'showcallinfo', 'no');
      tmpcfg = ft_setopt(tmpcfg, 'trackcallinfo', 'no');
      tmplayout = ft_prepare_layout(tmpcfg, data);
      layout.color = tmplayout.color;
    catch
      % it will fail if the data does not contain grad or elec
      ft_warning('cannot determine spatial colors');
    end
  end

elseif isequal(cfg.layout, 'vertical') || isequal(cfg.layout, 'horizontal')
  assert(iscell(cfg.channel), 'cfg.channel should be a cell-array of strings');
  nchan        = length(cfg.channel);
  layout.label = cfg.channel;
  
  % the width and height of the box are as specified
  % the distance between the channels is slightly larger
  switch cfg.layout
    case 'vertical'
      cfg.width  = ft_getopt(cfg, 'width',  0.8);
      cfg.height = ft_getopt(cfg, 'height', 0.8 * 1/(nchan+1+2));
    case 'horizontal'
      cfg.width  = ft_getopt(cfg, 'width',  0.8 * 1/(nchan+1+2));
      cfg.height = ft_getopt(cfg, 'height', 0.8);
  end
  
  for i=1:(nchan+2)
    switch cfg.layout
      case 'vertical'
        switch upper(cfg.direction)
          case 'TB'
            y = 1 - i*(cfg.height/0.8);
          case 'BT'
            y = 0 + i*(cfg.height/0.8);
          otherwise
            ft_error('invalid direction "%s" for "%s"', cfg.direction, cfg.layout);
        end
        x = 0.5*(cfg.width/0.8);
        layout.pos   (i,:) = [x y];
        layout.width (i,1) = cfg.width;
        layout.height(i,1) = cfg.height;
      case 'horizontal'
        switch upper(cfg.direction)
          case 'LR'
            x = 0 + i*(cfg.width/0.8);
          case 'RL'
            x = 1 - i*(cfg.width/0.8);
          otherwise
            ft_error('invalid direction "%s" for "%s"', cfg.direction, cfg.layout);
        end
        y = 0.5*(cfg.height/0.8);
        layout.pos   (i,:) = [x y];
        layout.width (i,1) = cfg.width;
        layout.height(i,1) = cfg.height;
    end
    if i==(nchan+1)
      layout.label{i}   = 'SCALE';
    elseif i==(nchan+2)
      layout.label{i}   = 'COMNT';
    end
  end
  
elseif isequal(cfg.layout, 'ordered')
  assert(iscell(cfg.channel), 'cfg.channel should be a valid set of channels');
  nchan        = length(cfg.channel);
  layout.label = cfg.channel;
  
  % the user can specify the number of columns and rows
  if isfield(cfg, 'columns') && ~isempty(cfg.columns)
    ncol = ft_getopt(cfg, 'columns');
  else
    ncol = nan; % wil be determined further down
  end
  if isfield(cfg, 'rows') && ~isempty(cfg.rows)
    nrow = ft_getopt(cfg, 'rows');
  else
    nrow = nan; % wil be determined further down
  end
  if isnan(ncol) && isnan(nrow)
    % the default is a more-or-less square arrangement
    ncol = ceil(sqrt(nchan))+1;
    nrow = ceil(sqrt(nchan))+1;
  elseif isnan(ncol)
    ncol = ceil(nchan/nrow);
  elseif isnan(nrow)
    nrow = ceil(nchan/ncol);
  end
  
  switch upper(cfg.direction)
    case 'LRTB'
      [Y, X] = ndgrid(1:nrow, 1:ncol);
      Y = flipud(Y);
      X = X';
      Y = Y';
    case 'RLTB'
      [Y, X] = ndgrid(1:nrow, 1:ncol);
      X = fliplr(X);
      Y = flipud(Y);
      X = X';
      Y = Y';
    case 'LRBT'
      [Y, X] = ndgrid(1:nrow, 1:ncol);
      X = X';
      Y = Y';
    case 'RLBT'
      [Y, X] = ndgrid(1:nrow, 1:ncol);
      X = fliplr(X);
      X = X';
      Y = Y';
    case 'TBLR'
      [Y, X] = ndgrid(1:nrow, 1:ncol);
      Y = flipud(Y);
    case 'TBRL'
      [Y, X] = ndgrid(1:nrow, 1:ncol);
      X = fliplr(X);
      Y = flipud(Y);
    case 'BTLR'
      [Y, X] = ndgrid(1:nrow, 1:ncol);
    case 'BTRL'
      [Y, X] = ndgrid(1:nrow, 1:ncol);
      X = fliplr(X);
    otherwise
      ft_error('invalid direction "%s" for "%s"', cfg.direction, cfg.layout);
  end
  
  cfg.width  = ft_getopt(cfg, 'width',  0.8 * 1/ncol);
  cfg.height = ft_getopt(cfg, 'height', 0.8 * 1/nrow);
  
  X = (X-1)*(cfg.width/0.8);
  Y = (Y-1)*(cfg.height/0.8);
  layout.pos = [X(:) Y(:)];
  layout.pos = layout.pos(1:nchan,:);
  
  layout.width  = ones(nchan,1) * cfg.width;
  layout.height = ones(nchan,1) * cfg.height;
  
  x = max(layout.pos(:,1));
  y = min(layout.pos(:,2)) - (cfg.height/0.8);
  scalepos = [x y];
  x = min(layout.pos(:,1));
  y = min(layout.pos(:,2)) - (cfg.height/0.8);
  comntpos = [x y];
  
  layout.label{end+1}  = 'SCALE';
  layout.pos(end+1,:)  = scalepos;
  layout.width(end+1)  = cfg.width;
  layout.height(end+1) = cfg.height;
  
  layout.label{end+1}  = 'COMNT';
  layout.pos(end+1,:)  = comntpos;
  layout.width(end+1)  = cfg.width;
  layout.height(end+1) = cfg.height;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % try to generate layout from other configuration options
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ischar(cfg.layout)
  [p, f, x] = fileparts(cfg.layout);
  
  if isempty(p) && isempty(x)
    % this is not a complete filename
    % check whether a corresponding lay or mat exists
    if exist([cfg.layout '.lay'], 'file')
      cfg.layout = [cfg.layout '.lay'];
      layout = ft_prepare_layout(cfg);
      return
    elseif exist([lower(cfg.layout) '.lay'], 'file')
      cfg.layout = [cfg.layout '.lay'];
      layout = ft_prepare_layout(cfg);
      return
    elseif exist([cfg.layout '.mat'], 'file')
      cfg.layout = [cfg.layout '.mat'];
      layout = ft_prepare_layout(cfg);
      return
    elseif exist([lower(cfg.layout) '.mat'], 'file')
      cfg.layout = [cfg.layout '.mat'];
      layout = ft_prepare_layout(cfg);
      return
    end
    
  elseif ft_filetype(cfg.layout, 'matlab')
    ft_info('reading layout from file %s\n', cfg.layout);
    if ~exist(cfg.layout, 'file')
      ft_error('the specified layout file %s was not found', cfg.layout);
    end
    layout = loadvar(cfg.layout,'layout');
    
  elseif ft_filetype(cfg.layout, 'layout')
    if exist(cfg.layout, 'file')
      ft_info('reading layout from file %s\n', cfg.layout);
      layout = readlay(cfg.layout);
    else
      [p, f] = fileparts(cfg.layout);
      ft_warning('the file "%s" was not found on your path, attempting "%s" instead', cfg.layout, fullfile(p, [f '.mat']));
      cfg.layout = fullfile(p, [f '.mat']);
      layout = ft_prepare_layout(cfg);
      return;
    end
    
  elseif ~ft_filetype(cfg.layout, 'layout')
    % assume that it points to an electrode file
    ft_info('creating layout from sensor description file %s\n', cfg.layout);
    sens = ft_read_sens(cfg.layout);
    layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
    
  end

elseif ~isempty(cfg.grad)
  if isstruct(cfg.grad)
    ft_info('creating layout from cfg.grad\n');
    sens = ft_datatype_sens(cfg.grad);
  else
    ft_info('creating layout from gradiometer file %s\n', cfg.grad);
    sens = ft_read_sens(cfg.grad, 'senstype', 'meg');
  end
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  
elseif isfield(data, 'grad') && isstruct(data.grad)
  ft_info('creating layout from data.grad\n');
  sens = ft_datatype_sens(data.grad);
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  
elseif ~isempty(cfg.elec)
  if isstruct(cfg.elec)
    ft_info('creating layout from cfg.elec\n');
    sens = ft_datatype_sens(cfg.elec);
  else
    ft_info('creating layout from electrode file %s\n', cfg.elec);
    sens = ft_read_sens(cfg.elec, 'senstype', 'eeg');
  end
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  
elseif isfield(data, 'elec') && isstruct(data.elec)
  ft_info('creating layout from data.elec\n');
  sens = ft_datatype_sens(data.elec);
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  
elseif ~isempty(cfg.opto)
  if isstruct(cfg.opto)
    ft_info('creating layout from cfg.opto\n');
    sens = ft_datatype_sens(cfg.opto);
  else
    ft_info('creating layout from optode file %s\n', cfg.opto);
    sens = ft_read_sens(cfg.opto, 'senstype', 'nirs');
  end
  if (hasdata)
    layout = opto2lay(sens, data.label, cfg.rotate, cfg.projection, cfg.viewpoint);
  else
    layout = opto2lay(sens, sens.label, cfg.rotate, cfg.projection, cfg.viewpoint);
  end
  
elseif isfield(data, 'opto') && isstruct(data.opto)
  ft_info('creating layout from data.opto\n');
  sens = ft_datatype_sens(data.opto);
  if (hasdata)
    layout = opto2lay(sens, data.label, cfg.rotate, cfg.projection, cfg.viewpoint);
  else
    layout = opto2lay(sens, sens.label, cfg.rotate, cfg.projection, cfg.viewpoint);
  end
  
elseif (~isempty(cfg.image) || ~isempty(cfg.mesh)) && isempty(cfg.layout)
  if ~isempty(cfg.image)
    % deal with image file
    ft_info('creating layout from image %s\n', cfg.image);
    [p, f, e] = fileparts(cfg.image);
    switch e
      case '.mat'
        img = loadvar(cfg.image);
      otherwise
        img = imread(cfg.image);
    end
    img = flip(img, 1); % in combination with "axis xy"
    
    figure
    
    if istrue(cfg.bw)
      % convert to greyscale image
      img = mean(img, 3);
      imagesc(img);
      colormap gray
    else
      % plot as RGB image
      image(img);
    end
    
  elseif ~isempty(cfg.mesh)
    if isfield(cfg.mesh, 'sulc')
      ft_plot_mesh(cfg.mesh, 'edgecolor', 'none', 'vertexcolor', cfg.mesh.sulc); colormap gray;
    else
      ft_plot_mesh(cfg.mesh, 'edgecolor', 'none');
    end
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
    'press the left mouse button to add another electrode\n' ...
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
          if ~isempty(cfg.image)
            h = image(img);
          else
            h = ft_plot_mesh(cfg.mesh, 'edgecolor', 'none', 'vertexcolor', cfg.mesh.sulc);
          end
          hold on
          axis equal
          axis off
          plot(pos(:,1), pos(:,2), 'b.');
          plot(pos(:,1), pos(:,2), 'yo');
        end
        
      case 'q'
        again = 0;
        
      otherwise
        ft_warning('invalid button (%d)', k);
    end
  end
  
  % get the interpolation mask
  polygon = {};
  thispolygon = 1;
  polygon{thispolygon} = zeros(0,2);
  maskhelp = [ ...
    '------------------------------------------------------------------------\n' ...
    'specify polygons for masking the topgraphic interpolation\n' ...
    'press the left mouse button to add another point to the current polygon\n' ...
    'press backspace on the keyboard to remove the last point\n' ...
    'press "c" on the keyboard to close this polygon and start with another\n' ...
    'press "q" on the keyboard to continue\n' ...
    ];
  again = true;
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
          if ~isempty(cfg.image)
            h = image(img);
          else
            h = ft_plot_mesh(cfg.mesh, 'edgecolor', 'none', 'vertexcolor', cfg.mesh.sulc);
          end
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
        ft_warning('invalid button (%d)', k);
    end
  end % while again
  % remember this set of polygons as the mask
  mask = polygon;
  
  % get the outline, e.g. head shape, nose, ears, sulci, etc.
  polygon = {};
  thispolygon = 1;
  polygon{thispolygon} = zeros(0,2);
  outlinehelp = [ ...
    '-----------------------------------------------------------------------------------\n' ...
    'specify polygons for adding outlines (e.g. head shape and sulci) to the layout\n' ...
    'press the left mouse button to add another point to the current polygon\n' ...
    'press backspace on the keyboard to remove the last point\n' ...
    'press "c" on the keyboard to close this polygon and start with another\n' ...
    'press "n" on the keyboard to start with another without closing the current polygon\n' ...
    'press "q" on the keyboard to continue\n' ...
    ];

  again = true;
  while again
    fprintf(outlinehelp);
    fprintf('\n');
    for i=1:length(polygon)
      fprintf('polygon %d has %d points\n', i, size(polygon{i},1));
    end
    
    try
      [x, y, k] = ginput(1);
    catch
      % this happens if the figure is closed
      return
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
          if ~isempty(cfg.image)
            image(img);
          else
            ft_plot_mesh(cfg.mesh, 'edgecolor', 'none', 'vertexcolor', cfg.mesh.sulc);
          end
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
        ft_warning('invalid button (%d)', k);
    end
  end % while again
  % remember this set of polygons as the outline
  outline = polygon;
  
  % convert the sensor positions into a layout structure
  layout.pos = pos;
  nchans = size(pos,1);
  for i=1:nchans
    layout.label{i,1} = num2str(i);
  end
  % compute the width and height for multiplotting
  d = dist(pos');
  for i=1:nchans
    d(i,i) = inf; % exclude the diagonal
  end
  mindist = min(d(:));
  layout.width  = ones(nchans,1) * mindist * 0.8;
  layout.height = ones(nchans,1) * mindist * 0.6;
  % add the polygons that describe the mask and outline
  layout.mask    = mask;
  layout.outline = outline;
  
  finalhelp = [ ...
    '-----------------------------------------------------------------------------------\n' ...
    'you should update the channel labels, and check the width and height in the output layout\n' ...
    ];
  fprintf(finalhelp);
  fprintf('\n');
  
else
  ft_error('no layout detected, please specify cfg.layout')
end

% make the subset as specified in cfg.channel
cfg.channel = ft_channelselection(cfg.channel, setdiff(layout.label, {'COMNT', 'SCALE'}, 'stable'));  % exclude COMNT and SCALE which are not really channels
chansel = match_str(layout.label, vertcat(cfg.channel(:), 'COMNT', 'SCALE'));                         % include COMNT and SCALE at the end, keep other channels in the order of the layout
% return the layout for the subset of channels
layout.pos    = layout.pos(chansel,:);
layout.label  = layout.label(chansel);
if strcmpi(cfg.style, '2d')
  % width and height only apply to the 2D layout
  layout.width  = layout.width(chansel);
  layout.height = layout.height(chansel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overrule the width and height when specified by the user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.width)
  layout.width(:) = cfg.width;
end
if ~isempty(cfg.height)
  layout.height(:) = cfg.height;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether the outline and mask are available, create them if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~isfield(layout, 'outline') || ~isfield(layout, 'mask')) && ~strcmpi(cfg.style, '3d')
  % the reason to check for style=3d rather than 2d is that cfg.style is also an option in ft_topoplotER and ft_topoplotTFR
  % the style option of that function easily "leaks" into here, causing the default 2d not to be selected at the top
  
  if strcmp(cfg.outline, 'circle') || strcmp(cfg.outline, 'doublecirclecross')
    % Scale the electrode positions to fit well within a unit circle
    ind_scale = find(strcmp('SCALE', layout.label));
    ind_comnt = find(strcmp('COMNT', layout.label));
    sel = setdiff(1:length(layout.label), [ind_scale ind_comnt]); % these are excluded for scaling
    x = layout.pos(sel,1);
    y = layout.pos(sel,2);
    if strcmp(cfg.center, 'yes')
      % the following centers all electrodes around zero
      xrange = range(x);
      yrange = range(y);
      shiftx = min(x);
      shifty = min(y);
    elseif strcmp(cfg.center, 'no')
      % the following prevent topography distortion in case electrodes are not evenly distributed over the whole head
      xrange = 2*( max(max(x),abs(min(x)) ));
      yrange = 2*( max(max(y),abs(min(y)) ));
      shiftx =   ( max(max(x),abs(min(x)) )).*sign(min(x));
      shifty =   ( max(max(y),abs(min(y)) )).*sign(min(y));
    end
    if xrange==0
      xrange = 1;
    end
    if yrange==0
      yrange = 1;
    end
    % First scale the width and height of the boxes for multiplotting
    layout.width  = layout.width./xrange;
    layout.height = layout.height./yrange;
    % Then shift and scale the electrode positions
    layout.pos(:,1) = 0.8*((layout.pos(:,1)-shiftx)/xrange-0.5);
    layout.pos(:,2) = 0.8*((layout.pos(:,2)-shifty)/yrange-0.5);
  end
  
  if ~isfield(layout, 'outline') && ischar(cfg.outline)
    switch cfg.outline
      case 'circle'
        layout.outline = outline_circle();
      case 'doublecirclecross'
        layout.outline = outline_doublecirclecross();
      case 'helmet'
        layout.outline = outline_helmet();
      case 'square'
        layout.outline = outline_square(layout);
      case 'convex'
        layout.outline = outline_convex(layout);
      case {'headshape', 'mri'}
        % the configuration should contain the headshape or mri
        % the (segmented) mri will be converted into a headshape on the fly
        hsoutline = outline_headshape(cfg, sens); % used for mask if possible
        layout.outline = hsoutline;
      otherwise
        layout.outline = {};
    end
  end
  
  if ~isfield(layout, 'mask') && ischar(cfg.mask)
    switch cfg.mask
      case 'circle'
        layout.mask = outline_circle();
        layout.mask = layout.mask(1); % the first is the outermost contour
      case 'helmet'
        layout.mask = outline_helmet();
        layout.mask = layout.mask(1); % the first is the outermost contour
      case 'extended'
        % take the outermost contour from the outline and extend it a bit
        layout.mask{1} = 1.2 * layout.outline{1}; 
      case 'square'
        layout.mask = outline_square(layout);
      case 'convex'
        layout.mask = outline_convex(layout);
      case {'headshape', 'mri'}
        % the configuration should contain the headshape or mri
        % the (segmented) mri will be converted into a headshape on the fly
        if isequal(cfg.mask, cfg.outline) && exist('hsoutline', 'var')
          layout.mask = hsoutline; % reuse the one from above
        else
          layout.mask = outline_headshape(cfg, sens);
        end
      otherwise
        layout.mask = {};
    end
  end
  
end % create outline and mask if style=2d


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the montage, e.g. convert from monopolar to bipolar channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(cfg.montage, 'no')
  Nold = length(cfg.montage.labelold);
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
if skipcomnt || ~isequal(cfg.commentpos, 'layout')
  % remove the comnt entry
  sel = find(strcmp('COMNT', layout.label));
  layout.label(sel)  = [];
  layout.pos(sel,:)  = [];
  layout.width(sel)  = [];
  layout.height(sel) = [];
end

if skipscale || ~isequal(cfg.scalepos, 'layout')
  % remove the scale entry
  sel = find(strcmp('SCALE', layout.label));
  layout.label(sel)  = [];
  layout.pos(sel,:)  = [];
  layout.width(sel)  = [];
  layout.height(sel) = [];
end

if (~skipcomnt || ~skipscale) && ~strcmpi(cfg.style, '3d')
  % this is used for the placement of the comment and scale
  pos = layout.pos;
  if isfield(layout, 'outline')
    pos = cat(1, pos, layout.outline{:});
  end
  if isfield(layout, 'mask')
    pos = cat(1, pos, layout.mask{:});
  end
  width  = mean(layout.width);
  height = mean(layout.height);
  middle = @(x) min(x) + (max(x)-min(x))/2;
end

if ~skipcomnt && ~any(strcmp('COMNT', layout.label)) && ~strcmpi(cfg.style, '3d') && ~isequal(cfg.commentpos, 'title')
  % add a placeholder for the comment in the desired location
  if strcmp(cfg.commentpos, 'layout')
    cfg.commentpos = 'leftbottom'; % set the default position
  end
  if strcmp(cfg.commentpos, 'lefttop')
    layout.pos(end+1,:) = [min(pos(:,1))-width/2 max(pos(:,2))+height/2];
  elseif strcmp(cfg.commentpos, 'leftbottom')
    layout.pos(end+1,:) = [min(pos(:,1))-width/2 min(pos(:,2))-height/2];
  elseif strcmp(cfg.commentpos, 'middletop')
    layout.pos(end+1,:) = [middle(pos(:,1)) max(pos(:,2))+height/2];
  elseif strcmp(cfg.commentpos, 'middlebottom')
    layout.pos(end+1,:) = [middle(pos(:,1)) min(pos(:,2))-height/2];
  elseif strcmp(cfg.commentpos, 'righttop')
    layout.pos(end+1,:) = [max(pos(:,1))+width/2 max(pos(:,2))+height/2];
  elseif strcmp(cfg.commentpos, 'rightbottom')
    layout.pos(end+1,:) = [max(pos(:,1))+width/2 min(pos(:,2))-height/2];
  elseif isnumeric(cfg.commentpos)
    layout.pos(end+1,:) = cfg.commentpos;
  else
    ft_error('invalid specification of cfg.commentpos');
  end
  layout.label{end+1}  = 'COMNT';
  layout.width(end+1)  = width;
  layout.height(end+1) = height;
end

if ~skipscale && ~any(strcmp('SCALE', layout.label)) && ~strcmpi(cfg.style, '3d')
  % add a placeholder for the scale in the desired location
  if strcmp(cfg.scalepos, 'layout')
    cfg.scalepos = 'rightbottom'; % set the default position
  end
  if strcmp(cfg.scalepos, 'lefttop')
    layout.pos(end+1,:) = [min(pos(:,1))-width/2 max(pos(:,2))+height/2];
  elseif strcmp(cfg.scalepos, 'leftbottom')
    layout.pos(end+1,:) = [min(pos(:,1))-width/2 min(pos(:,2))-height/2];
  elseif strcmp(cfg.scalepos, 'middletop')
    layout.pos(end+1,:) = [middle(pos(:,1)) max(pos(:,2))+height/2];
  elseif strcmp(cfg.scalepos, 'middlebottom')
    layout.pos(end+1,:) = [middle(pos(:,1)) min(pos(:,2))-height/2];
  elseif strcmp(cfg.scalepos, 'righttop')
    layout.pos(end+1,:) = [max(pos(:,1))+width/2 max(pos(:,2))+height/2];
  elseif strcmp(cfg.scalepos, 'rightbottom')
    layout.pos(end+1,:) = [max(pos(:,1))+width/2 min(pos(:,2))-height/2];
  elseif isnumeric(cfg.scalepos)
    layout.pos(end+1,:) = cfg.scalepos;
  else
    ft_error('invalid specification of cfg.scalepos');
  end
  layout.label{end+1}  = 'SCALE';
  layout.width(end+1)  = width;
  layout.height(end+1) = height;
end

% these should be represented in a column vector (see bug 1909 -roevdmei)
layout.label = layout.label(:);
% the width and height are not present in a 3D layout as used in SPM
if ~strcmpi(cfg.style, '3d')
  layout.width  = layout.width(:);
  layout.height = layout.height(:);
end

% to plot the layout for debugging, you can use this code snippet
if strcmp(cfg.feedback, 'yes')
  if strcmpi(cfg.style, '3d')
    ft_error('graphical feedback is not implemented for a 3d layout');
  else
    ft_plot_layout(layout);
  end
end

% to write the layout to a .mat or text file, you can use this code snippet
if ~isempty(cfg.output) && ~strcmpi(cfg.style, '3d')
  ft_info('writing layout to ''%s''\n', cfg.output);
  if strcmpi(cfg.output((end-3):end), '.mat')
    save(cfg.output,'layout');
  else
    fid = fopen(cfg.output, 'wt');
    for i=1:numel(layout.label)
      fprintf(fid, '%d %f %f %f %f %s\n', i, layout.pos(i,1), layout.pos(i,2), ...
        layout.width(i), layout.height(i), layout.label{i});
    end
    fclose(fid);
  end
elseif ~isempty(cfg.output) && strcmpi(cfg.style, '3d')
  % the layout file format does not support 3D positions, furthermore for
  % a 3D layout the width and height are currently set to NaN
  ft_error('writing a 3D layout to an output file is not supported');
end

if isequal(cfg.color, 'spatial') && ~isfield(layout, 'color')
  % create a channel-specific RGB value based on their X/Y positions and a
  % computed Z position, assuming the channels on the positive half sphere
  sel = match_str(layout.label, setdiff(layout.label, {'COMNT';'SCALE'}));
  xy  = sqrt(sum(layout.pos(sel,:).^2,2));
  z   = sqrt(max(xy).^2 - xy.^2);
  xyz = [layout.pos(sel,:) z];
  xyz = xyz - min(xyz, [], 1);
  rgb = xyz./(max(xyz, [], 1)+10*eps);
  layout.color = ones(numel(layout.label), 3);
  layout.color(sel, :) = rgb;
elseif size(cfg.color,1) == numel(layout.label)
  % copy the specified colors over
  layout.color = cfg.color;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble history layout
ft_postamble savevar layout


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% read the layout information from the ascii file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = readlay(filename)
if ~exist(filename, 'file')
  ft_error('could not open layout file "%s"', filename);
end
fid=fopen(filename);
lay_string=fread(fid,inf,'char=>char')';
fclose(fid);

% pattern to match is integer, than 4 numeric values followed by a
% string that can contain whitespaces and plus characters, followed by
% newline
integer='(\d+)';
float='([\d\.-]+)';
space='\s+';
channel_label='([\w \t\r\f\v\-\+]+)';
single_newline='\n';

pat=[integer, space, ...
  float, space, ...
  float, space, ...
  float, space, ...
  float, space, ...
  channel_label, single_newline];

matches=regexp(sprintf('%s\n',lay_string),pat,'tokens');

% convert to (nchannel x 6) matrix
layout_matrix=cat(1,matches{:});

% convert values in first five columns to numeric
num_values_cell=layout_matrix(:,1:5)';
str_values=sprintf('%s %s %s %s %s; ', num_values_cell{:});
num_values=str2num(str_values);

% store layout information (omit channel number in first column)
layout.pos    = num_values(:,2:3);
layout.width  = num_values(:,4);
layout.height = num_values(:,5);

% trim whitespace around channel names
label=layout_matrix(:,6);
label=regexprep(label,'^\s*','');
label=regexprep(label,'\s*$','');
layout.label  = label;
return % function readlay


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% convert 3D electrode positions into 2D layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = sens2lay(sens, rotatez, projmethod, style, overlap, viewpoint, boxchannel)

% remove the balancing from the sensor definition, e.g. 3rd order gradients, PCA-cleaned data or ICA projections
% this not only removed the linear projections, but also ensures that the channel labels are correctly named

if isfield(sens, 'chanposold')
  chanposold = sens.chanposold;
else
  chanposold = [];
end
if isfield(sens, 'balance') && ~strcmp(sens.balance.current, 'none')
  sens = undobalancing(sens);
  if size(chanposold, 1) == numel(sens.label)
    sens.chanpos = chanposold;
  end
  % In case not all the locations have NaNs it might still be useful to plot them
  % But perhaps it'd be better to have any
elseif any(all(isnan(sens.chanpos)))
  [sel1, sel2] = match_str(sens.label, sens.labelold);
  sens.chanpos = chanposold(sel2, :);
  sens.label   = sens.labelold(sel2);
end

ft_info('creating layout for %s system\n', ft_senstype(sens));

% apply rotation, but only if viewpoint is not used specifically
if isempty(viewpoint) && isempty(rotatez)
  if isfield(sens, 'coordsys')
    % the x-axis to the right of the screen and the y-axis is to the top of the screen
    % the nose is usually plotted toward the top of the screen
    sens = ft_convert_coordsys(sens, 'ras', 0);
  elseif isfield(sens, 'type') && startsWith(sens.type, 'ctf')
    sens.coordsys = 'ctf';
    sens = ft_convert_coordsys(sens, 'ras', 0);
  elseif isfield(sens, 'type') && startsWith(sens.type, 'bti')
    sens.coordsys = 'bti';
    sens = ft_convert_coordsys(sens, 'ras', 0);
  elseif isfield(sens, 'type') && startsWith(sens.type, '4d')
    sens.coordsys = '4d';
    sens = ft_convert_coordsys(sens, 'ras', 0);
  elseif isfield(sens, 'type') && startsWith(sens.type, 'yokogawa')
    sens.coordsys = 'yokogawa';
    sens = ft_convert_coordsys(sens, 'ras', 0);
  elseif isfield(sens, 'type') && startsWith(sens.type, 'neuromag')
    sens.coordsys = 'neuromag';
    sens = ft_convert_coordsys(sens, 'ras', 0);
  elseif isfield(sens, 'type') && startsWith(sens.type, 'itab')
    sens.coordsys = 'itab';
    sens = ft_convert_coordsys(sens, 'ras', 0);
  else
    ft_warning('use cfg.rotate to specify the correct rotation of the sensors');
  end
end

% determine the 3D channel positions
pos   = sens.chanpos;
label = sens.label;

if ~isempty(rotatez)
  % apply the rotation around the z-axis
  pos = ft_warp_apply(rotate([0 0 rotatez]), pos, 'homogenous');
end

if strcmpi(style, '3d')
  layout.pos   = pos;
  layout.label = label;
  
else
  if isempty(viewpoint)
    % projection other than viewpoint-specific orthographic projection is requested, use elproj
    prj = elproj(pos, projmethod);
    
  else
    % apply viewpoint-specific orthographic projection
    
    % determine auto view if requested
    if strcmp(viewpoint, 'auto')
      % simple automatic determination of the 'ideal' viewpoint
      % first, depth or not: if Xvar (l/r axis) is bigger than both Yvar (post/ant axis) and Zvar (top/bottom axis), it's a depth
      % if yes, superior (screw inferior) is more appriorate if Yvar > Zvar, otherwise posterior (screw anterior)
      % if no, it's left/right, sign of mean(X) indicates which side the grid is on (note, for interhemispheric grids, both left/right (doenst) work)
      posvar = var(pos);
      if (posvar(1)>posvar(2)) && (posvar(1)>posvar(3)) % if they're roughly equal, it's likely a diagonal depth, and any view would (not) work
        if posvar(2)>posvar(3)
          viewpoint = 'superior';
        else
          viewpoint = 'posterior';
        end
      else
        if sign(mean(pos(:,1))) == -1
          viewpoint = 'left';
        else
          viewpoint = 'right';
        end
      end
    end
    
    % project 3D points to 2D
    prj = getorthoviewpos(pos, sens.coordsys, viewpoint);
  end % if viewpoint
  
  % this copy will be used to determine the minimum distance between channels
  % we need a copy because prj retains the original positions, and
  % prjForDist might need to be changed if the user wants to keep
  % overlapping channels
  % also subselect channels for computing width/height if requested by boxchannel
  boxchannel = ft_channelselection(boxchannel, label);
  boxchansel = match_str(label, boxchannel);
  prjForDist = prj(boxchansel,:);
  
  % check whether many channels occupy identical positions, if so shift them around if requested
  if size(unique(prjForDist,'rows'),1) / size(prjForDist,1) < 0.8
    if strcmp(overlap, 'shift')
      ft_warning('the specified sensor configuration has many overlapping channels, creating a layout by shifting them around - use a template layout for better control over the positioning');
      prj = shiftxy(prj', 0.2)';
      prjForDist = prj(boxchansel,:);
    elseif strcmp(overlap, 'no')
      ft_error('the specified sensor configuration has many overlapping channels, you specified not to allow that');
    elseif strcmp(overlap, 'keep')
      prjForDist = unique(prj(boxchansel,:), 'rows');
    else
      ft_error('unknown value for cfg.overlap = ''%s''', overlap);
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
  % also exists for some .lay files.
  if any(strmatch('Fid', label(boxchansel)))
    tmpsel = strmatch('Fid', label(boxchansel));
    d(tmpsel, :) = inf;
    d(:, tmpsel) = inf;
  end
  
  if any(isfinite(d(:)))
    % take mindist as the median of the first quartile of closest channel pairs with non-zero distance
    mindist = min(d); % get closest neighbour for all channels
    mindist = sort(mindist(mindist>1e-6),'ascend');
    mindist = mindist(1:round(numel(label(boxchansel))/4));
    mindist = median(mindist);
    %%% /workaround - make a safe guess to detect iEEG until a better solution is found
    if any(strcmp(ft_senstype(sens),{'eeg','unknown'})) && ~isempty(viewpoint) && (numel(boxchannel) ~= numel(label))
      mindist = min(d);
      mindist = mindist(mindist>1e-6); % allows for substantially more overlap than the above
      mindist = median(mindist);
    end
    %%% \workaround - make a safe guess to detect iEEG until a better solution is found
  else
    mindist = eps; % not sure this is a good value, but it's just to prevent crashes when the EEG sensor definition is meaningless
  end
  
  layout.pos    = prj;
  layout.label  = label;
  layout.width  = ones(numel(label),1) * mindist * 0.8;
  layout.height = ones(numel(label),1) * mindist * 0.6;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% convert 2D optode positions into 2D layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = opto2lay(opto, label, rotatez, projmethod, viewpoint)

if isempty(rotatez)
  rotatez = 90;
end
  
% NIRS channels are named as 'RxY-TxZ [wavelength]' or as 'RxY-TxZ [chromophore]'
[rxnames, rem] = strtok(label, {'-', ' '});
[txnames, rem] = strtok(rem,   {'-', ' '});

pos = nan(numel(label),3);
for i=1:numel(label)
  if isfield(opto, 'chanpos') && ismember(label{i}, opto.label)
    % there is an exact match
    chanid = strcmp(opto.label, label{i});
    pos(i,:) = opto.chanpos(chanid,:);
  elseif isfield(opto, 'chanpos') && any(startsWith(opto.label, [rxnames{i} '-' txnames{i}]))
    % the first part matches with 'RxY-TxZ'
    chanid = startsWith(opto.label, [rxnames{i} '-' txnames{i}]);
    pos(i,:) = mean(opto.chanpos(chanid,:), 1); % there will usually be two matches
  elseif isfield(opto, 'chanpos') && any(startsWith(opto.label, [txnames{i} '-' rxnames{i}]))
    % the first part matches with 'TxY-RxZ'
    chanid = startsWith(opto.label, [txnames{i} '-' rxnames{i}]);
    pos(i,:) = mean(opto.chanpos(chanid,:), 1); % there will usually be two matches
  elseif ismember(rxnames(i), opto.optolabel) && ismember(txnames{i}, opto.optolabel)
    % create positions halfway between the transmitter and receiver
    rxid = strcmp(opto.optolabel, rxnames(i));
    txid = strcmp(opto.optolabel, txnames(i));
    pos(i, :) = opto.optopos(rxid, :)/2 + opto.optopos(txid, :)/2;
  end
end

% remove the channels without position, like AUX
sel = ~any(isnan(pos), 2);
pos = pos(sel,:);
label = label(sel);

% apply the rotation around the z-axis
pos = ft_warp_apply(rotate([0 0 rotatez]), pos, 'homogenous');

% project 3D points onto 2D plane
if all(pos(:,3)==0)
  ft_notice('not applying 2D projection');
  pos = pos(:,1:2);
elseif isempty(viewpoint)
  pos = elproj(pos, projmethod);
else
  pos = getorthoviewpos(pos, opto.coordsys, viewpoint);
end

% compute the distances between all channel pairs
dist = zeros(numel(label));
for i=1:numel(label)
  for j=1:numel(label)
    dist(i,j) = norm(pos(i,:)-pos(j,:));
  end
end
dist(dist==0) = inf; % ignore all zeros
mindist = median(min(dist,[], 2))/2;

if isinf(mindist) || isnan(mindist)
  % this is needed e.g. in case there is only one channel
  mindist = 1;
end

% note that the width and height can be overruled elsewhere with cfg.width and cfg.height
ft_notice('estimated channel width and height is %.4f', mindist);

% start with an empty layout
layout = [];
layout.pos    = pos;
layout.label  = label;
layout.width  = mindist*ones(numel(label),1);
layout.height = mindist*ones(numel(label),1);

% prevent the circle-with-ears-and-nose to be added
layout.outline = {};

% construct a mask for topographic interpolation
pos1 = layout.pos; pos1(:,1) = pos1(:,1) - layout.width; pos1(:,2) = pos1(:,2) - layout.height;
pos2 = layout.pos; pos2(:,1) = pos2(:,1) - layout.width; pos2(:,2) = pos2(:,2) + layout.height;
pos3 = layout.pos; pos3(:,1) = pos3(:,1) + layout.width; pos3(:,2) = pos3(:,2) - layout.height;
pos4 = layout.pos; pos4(:,1) = pos4(:,1) + layout.width; pos4(:,2) = pos4(:,2) + layout.height;
pos = [pos1; pos2; pos3; pos4];

indx = convhull(pos);
layout.mask{1} = pos(indx,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% shift 2D positions around so that the minimum distance between any pair
% is mindist
%
% Credit for this code goes to Laurence Hunt at UCL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xy = shiftxy(xy, mindist)

x = xy(1,:);
y = xy(2,:);

l=1;
i=1; % filler
mindist = mindist/0.999; % limits the number of loops
while (~isempty(i) && l<50)
  xdiff = repmat(x,length(x),1) - repmat(x',1,length(x));
  ydiff = repmat(y,length(y),1) - repmat(y',1,length(y));
  xydist= sqrt(xdiff.^2 + ydiff.^2); %euclidean distance between all sensor pairs
  
  [i,j] = find(xydist<mindist*0.999);
  rm=(i<=j); i(rm)=[]; j(rm)=[]; %only look at i>j
  
  for m = 1:length(i)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate an outline consisting of a unit-diameter circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_circle()
% create the default "circle with triangle" to resemble the head
% note that the electrode positions should be scaled accordingly
rmax  = 0.5;
l     = 0:2*pi/100:2*pi;
HeadX = cos(l).*rmax;
HeadY = sin(l).*rmax;
NoseX = [0.18*rmax 0 -0.18*rmax];
NoseY = [rmax-.004 rmax*1.15 rmax-.004];
EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
% Define the outline of the head, ears and nose
outline{1} = [HeadX(:) HeadY(:)];
outline{2} = [NoseX(:) NoseY(:)];
outline{3} = [ EarX(:)  EarY(:)];
outline{4} = [-EarX(:)  EarY(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate an outline consisting of a unit-diameter circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_doublecirclecross()
% create a double circle with a cross over the vertex
% note that the electrode positions should be scaled accordingly
rmax  = 0.5;
l     = 0:2*pi/100:2*pi;
HeadX = cos(l).*rmax;
HeadY = sin(l).*rmax;
NoseX = [0.18*rmax 0 -0.18*rmax];
NoseY = [rmax-.004 rmax*1.15 rmax-.004];
% Define the outline of the head with the cross
outline{1} = [HeadX(:) HeadY(:)]*1.0; % outer circle
outline{2} = [HeadX(:) HeadY(:)]*0.8; % inner circle
outline{3} = [NoseX(:) NoseY(:)];     % nose
outline{4} = [-0.5 0 ; +0.5 0];       % horizontal line
outline{5} = [0 -0.5 ; 0 +0.5];       % vertical line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate an outline consisting of default helmet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_helmet()
% this is taken from the layout for the FieldLine helmet, which in turn is taken from the CTF helmet
outline{1} = 5.5*[0.00331272465437779 0.470471964285714;-0.0120810253456222 0.470149464285714;-0.0274522753456222 0.469208214285714;-0.0582135253456222 0.465743214285714;-0.0812047753456221 0.461854464285714;-0.0831247753456221 0.461430714285714;-0.0850110253456223 0.460924464285714;-0.0859335253456222 0.460628214285714;-0.0868447753456222 0.460301964285714;-0.0877372753456222 0.459938214285714;-0.116811025345622 0.445429464285714;-0.131826025345622 0.437190714285714;-0.158338525345622 0.420795714285714;-0.161852275345622 0.418335714285714;-0.189647275345622 0.394980714285714;-0.212271025345622 0.373796964285714;-0.213921025345622 0.371974464285714;-0.214367275345622 0.371531964285714;-0.214588525345622 0.371329464285714;-0.214817275345622 0.371138214285714;-0.215046025345622 0.370965714285714;-0.215271025345622 0.370815714285714;-0.215383525345622 0.370748214285714;-0.215499775345622 0.370688214285714;-0.215627275345622 0.370628214285714;-0.215758525345622 0.370575714285714;-0.215889775345622 0.370534464285714;-0.216021025345622 0.370500714285714;-0.216298525345622 0.370455714285714;-0.216437275345622 0.370444464285714;-0.216579775345622 0.370436964285714;-0.216722275345622 0.370440714285714;-0.216864775345622 0.370448214285714;-0.217011025345622 0.370459464285714;-0.217303525345622 0.370500714285714;-0.217453525345622 0.370526964285714;-0.218057275345622 0.370673214285714;-0.265164775345622 0.384394464285714;-0.265629775345622 0.384488214285714;-0.266091025345622 0.384551964285714;-0.266544775345622 0.384589464285714;-0.266994775345622 0.384604464285714;-0.267441025345622 0.384596964285714;-0.268333525345622 0.384533214285714;-0.269233525345622 0.384405714285714;-0.307813525345622 0.377423214285714;-0.308031025345622 0.377363214285714;-0.308244775345622 0.377295714285714;-0.308451025345622 0.377216964285714;-0.308646025345622 0.377130714285714;-0.308837275345622 0.377036964285714;-0.309021025345622 0.376935714285714;-0.309197275345622 0.376819464285714;-0.309534775345622 0.376571964285714;-0.309853525345622 0.376298214285714;-0.310157275345622 0.376001964285714;-0.310303525345622 0.375844464285714;-0.310581025345622 0.375518214285714;-0.311113525345622 0.374820714285714;-0.312144775345622 0.373343214285714;-0.342174775345622 0.336293214285714;-0.356443525345622 0.316718214285714;-0.374994775345622 0.288619464285714;-0.392019775345622 0.259568214285714;-0.407477275345622 0.229654464285714;-0.421326025345622 0.198911964285714;-0.433442275345622 0.167404464285714;-0.443769775345622 0.135259464285714;-0.452244775345622 0.102585714285714;-0.459924775345622 0.0629069642857141;-0.464323525345622 0.0290219642857141;-0.465594775345622 0.012424464285714;-0.465681025345622 -0.00410928571428595;-0.464867275345622 -0.0206917857142859;-0.463476025345622 -0.037146785714286;-0.459418525345622 -0.069794285714286;-0.456302275345622 -0.087348035714286;-0.452556025345622 -0.104770535714286;-0.448172275345622 -0.122031785714286;-0.443169775345622 -0.139120535714286;-0.437548525345622 -0.156010535714286;-0.431312275345622 -0.172675535714286;-0.424483525345622 -0.189108035714286;-0.417054775345622 -0.205281785714286;-0.409041025345622 -0.221174285714286;-0.400449775345622 -0.236763035714286;-0.391288525345622 -0.252033035714286;-0.381561025345622 -0.266961785714286;-0.371278525345622 -0.281530535714286;-0.360444775345622 -0.295716785714286;-0.349074775345622 -0.309501785714286;-0.337194775345622 -0.322836785714286;-0.324827275345622 -0.335710535714286;-0.311991025345622 -0.348100535714286;-0.298712275345622 -0.359991785714286;-0.284994775345622 -0.371373035714286;-0.270864775345622 -0.382233035714286;-0.256344775345622 -0.392545535714286;-0.241453525345622 -0.402310535714286;-0.226194775345622 -0.411501785714286;-0.210602275345622 -0.420108035714286;-0.194687275345622 -0.428114285714286;-0.178472275345622 -0.435505535714286;-0.161968525345622 -0.442270535714286;-0.145206025345622 -0.448390535714286;-0.128184775345622 -0.453854285714286;-0.112352275345622 -0.458279285714286;-0.0830835253456223 -0.464856785714286;-0.0534360253456222 -0.469495535714286;-0.0235410253456222 -0.472191785714286;0.00647397465437775 -0.472960535714286;0.0364739746543777 -0.471783035714286;0.0663277246543778 -0.468674285714286;0.0959039746543778 -0.463630535714286;0.110543974654378 -0.460386785714286;0.126466474654378 -0.456265535714286;0.143588974654378 -0.451131785714286;0.160478974654378 -0.445330535714286;0.177113974654378 -0.438880535714286;0.193475224654378 -0.431789285714286;0.209547724654378 -0.424075535714286;0.225308974654378 -0.415761785714286;0.240743974654378 -0.406844285714286;0.255826474654378 -0.397353035714286;0.270545224654378 -0.387299285714286;0.284873974654378 -0.376698035714286;0.298801474654378 -0.365556785714286;0.312305224654378 -0.353898035714286;0.325366474654378 -0.341733035714286;0.337970224654378 -0.329080535714286;0.350090224654378 -0.315948035714286;0.361715224654378 -0.302354285714286;0.372811474654378 -0.288336785714286;0.383363974654378 -0.273929285714286;0.393368974654378 -0.259146785714286;0.402811474654378 -0.244015535714286;0.411687724654378 -0.228550535714286;0.419990224654378 -0.212781785714286;0.427707724654378 -0.196713035714286;0.434832724654378 -0.180378035714286;0.441361474654378 -0.163791785714286;0.447278974654378 -0.146980535714286;0.452585224654378 -0.129951785714286;0.457265224654378 -0.112739285714286;0.461315224654378 -0.0953542857142859;0.464723974654378 -0.0778230357142859;0.468863974654378 -0.0491730357142859;0.471601474654378 -0.0151492857142859;0.472032724654378 0.00191321428571403;0.471773974654378 0.0188857142857141;0.470588974654378 0.035306964285714;0.468515224654378 0.051728214285714;0.462372724654378 0.086348214285714;0.454790224654378 0.119269464285714;0.445351474654378 0.151706964285714;0.434105224654378 0.183563214285714;0.421096474654378 0.214721964285714;0.406433974654378 0.245018214285714;0.390158974654378 0.274478214285714;0.372342724654378 0.303038214285714;0.344903974654378 0.341509464285714;0.325103974654378 0.366308214285714;0.316085224654378 0.376808214285714;0.315886474654378 0.376939464285714;0.315672724654378 0.377063214285714;0.315455224654378 0.377179464285714;0.314986474654378 0.377385714285714;0.314742724654378 0.377479464285714;0.314228974654378 0.377644464285714;0.313688974654378 0.377790714285714;0.312552724654378 0.378023214285714;0.308015224654378 0.378596964285714;0.276117724654378 0.384480714285714;0.275056474654378 0.384630714285714;0.274535224654378 0.384679464285714;0.274017724654378 0.384701964285714;0.273511474654378 0.384701964285714;0.273260224654378 0.384694464285714;0.272765224654378 0.384649464285714;0.272363974654378 0.384593214285714;0.271745224654378 0.384461964285714;0.271122724654378 0.384304464285714;0.227900224654378 0.371670714285714;0.225391474654378 0.370733214285714;0.224716474654378 0.370511964285714;0.224390224654378 0.370421964285714;0.224071474654378 0.370350714285714;0.223767724654378 0.370298214285714;0.223621474654378 0.370279464285714;0.223478974654378 0.370268214285714;0.223340224654378 0.370264464285714;0.223208974654378 0.370264464285714;0.223081474654378 0.370275714285714;0.222968974654378 0.370290714285714;0.222852724654378 0.370320714285714;0.222740224654378 0.370354464285714;0.222623974654378 0.370399464285714;0.222507724654378 0.370451964285714;0.222391474654378 0.370511964285714;0.222155224654378 0.370654464285714;0.221922724654378 0.370823214285714;0.221686474654378 0.371010714285714;0.221221474654378 0.371434464285714;0.219541474654378 0.373193214285714;0.198893974654378 0.392655714285714;0.168526474654378 0.418298214285714;0.164982724654378 0.420788214285714;0.144886474654378 0.433380714285714;0.129991474654378 0.441968214285714;0.0952664746543777 0.459514464285714;0.0934402246543777 0.460238214285714;0.0915877246543778 0.460838214285714;0.0897052246543778 0.461340714285714;0.0743714746543778 0.464258214285714;0.0523402246543778 0.467351964285714;0.0215677246543778 0.470044464285714;0.00331272465437779 0.470471964285714]; 
outline{2} = 5.5*[4.5166850393727e-05 0.453549861023622;0.00402785185039376 0.453413916023622;0.00800550185039377 0.453190256023622;0.0119725518503937 0.452875701023622;0.0159268818503938 0.452475551023622;0.0198708768503937 0.451986626023622;0.0237992368503938 0.451413166023622;0.0277172618503938 0.450753316023622;0.0316249518503938 0.450008136023622;0.0355172718503938 0.449181071023622;0.0393989918503938 0.448270001023622;0.0432653418503938 0.447277576023622;0.0471160568503938 0.446203266023622;0.0509540518503938 0.445047866023622;0.0547764118503937 0.443813231023622;0.0585857868503937 0.442499891023622;0.0623795268503937 0.441108641023622;0.0661578968503938 0.439640011023622;0.0699206318503938 0.438094796023622;0.0736677318503938 0.436472996023622;0.0773941618503938 0.434778056023622;0.0811078718503938 0.433007591023622;0.0848032968503938 0.431166106023622;0.0884807018503938 0.429249891023622;0.0921400868503938 0.427262126023622;0.0957811868503937 0.425204401023622;0.0994042668503938 0.423075656023622;0.103009061850394 0.420878276023622;0.106593451850394 0.418612526023622;0.110156906850394 0.416279201023622;0.113702341850394 0.413878831023622;0.117227106850394 0.411412211023622;0.120730936850394 0.408880931023622;0.124214361850394 0.406284991023622;0.127676851850394 0.403624391023622;0.131118936850394 0.400903106023622;0.134537436850394 0.398119546023622;0.137935266850394 0.395272386023622;0.141307391850394 0.392367456023622;0.144658581850394 0.389401046023622;0.147986451850394 0.386376866023622;0.151288616850394 0.383294386023622;0.154570111850394 0.380156256023622;0.157825636850394 0.376960091023622;0.161055191850394 0.373709071023622;0.164264341850394 0.370402931023622;0.167442221850394 0.367040876023622;0.170599431850394 0.363628471023622;0.173730671850394 0.360162801023622;0.176831171850394 0.356646781023622;0.179908086850394 0.353079086023622;0.182959296850394 0.349460776023622;0.185981886850394 0.345794766023622;0.188976121850394 0.342079201023622;0.191944386850394 0.338316201023622;0.194881646850394 0.334507886023622;0.197792936850394 0.330652666023622;0.200673221850394 0.326752396023622;0.203525151850394 0.322808136023622;0.206348461850394 0.318819886023622;0.209138381850394 0.314788706023622;0.211902066850394 0.310717246023622;0.214635011850394 0.306603651023622;0.217334036850394 0.302449511023622;0.220002056850394 0.298255621023622;0.222641721850394 0.294023836023622;0.225245081850394 0.289754686023622;0.227820086850394 0.285448171023622;0.230361436850394 0.281105086023622;0.232868866850394 0.276727551023622;0.235342906850394 0.272314506023622;0.237785941850394 0.267866746023622;0.240192671850394 0.263387716023622;0.242563096850394 0.258875561023622;0.244902781850394 0.254332136023622;0.247205896850394 0.249757441023622;0.249472971850394 0.245153066023622;0.251706126850394 0.240521661023622;0.253903241850394 0.235860311023622;0.256061666850394 0.231171931023622;0.258186436850394 0.226456256023622;0.260274636850394 0.221714876023622;0.262324411850394 0.216951501023622;0.264335231850394 0.212161891023622;0.266309746850394 0.207347901023622;0.268245571850394 0.202513241023622;0.270139791850394 0.197655261023622;0.271997971850394 0.192776876023622;0.273817196850394 0.187878616023622;0.275595081850394 0.182959951023622;0.277334011850394 0.178024326023622;0.279034251850394 0.173070416023622;0.280693151850394 0.168099811023622;0.282310711850394 0.163112511023622;0.283886666850394 0.158109576023622;0.285418631850394 0.153093656023622;0.286911906850394 0.148062631023622;0.288363576850394 0.143019681023622;0.289771521850394 0.137964806023622;0.291132561850394 0.132898536023622;0.292454911850394 0.127821931023622;0.293733271850394 0.122733666023622;0.294967641850394 0.117638776023622;0.296158286850394 0.112534876023622;0.297302026850394 0.107424616023622;0.298401776850394 0.102306671023622;0.299455151850394 0.0971821010236217;0.300464536850394 0.0920548810236217;0.301427016850394 0.0869226260236217;0.302343121850394 0.0817877210236217;0.303212851850394 0.0766501660236217;0.304035676850394 0.0715112860236218;0.304812126850394 0.0663713460236217;0.305541936850394 0.0612311410236217;0.306219806850394 0.0560935860236217;0.306851301850394 0.0509565610236217;0.307436156850394 0.0458216560236218;0.307969336850394 0.0406891360236217;0.308455876850394 0.0355608560236217;0.308890476850394 0.0304376110236217;0.309276316850394 0.0253194010236217;0.309609951850394 0.0202104660236217;0.309894561850394 0.0151052410236217;0.310127496850394 0.0100106160236218;0.310311406850394 0.00492341102362175;0.310443376850394 -0.000153193976378283;0.310521021850394 -0.00521946397637828;0.310546726850394 -0.0102743389763783;0.310521021850394 -0.0153029789763782;0.310443376850394 -0.0202905439763783;0.310311406850394 -0.0252338539763782;0.310127496850394 -0.0301355589763783;0.309894561850394 -0.0349961889763783;0.309609951850394 -0.0398178639763782;0.309276316850394 -0.0445931639763783;0.308890476850394 -0.0493295089763783;0.308455876850394 -0.0540215989763783;0.307969336850394 -0.0586776489763783;0.307436156850394 -0.0632873239763783;0.306851301850394 -0.0678580439763783;0.306219806850394 -0.0723845089763782;0.305541936850394 -0.0768698989763782;0.304812126850394 -0.0813136839763783;0.304035676850394 -0.0857190439763783;0.303212851850394 -0.0900822689763783;0.302343121850394 -0.0944020339763783;0.301427016850394 -0.0986804589763782;0.300464536850394 -0.102922313976378;0.299455151850394 -0.107120178976378;0.298401776850394 -0.111276703976378;0.297302026850394 -0.115389503976378;0.296158286850394 -0.119465468976378;0.294967641850394 -0.123495588976378;0.293733271850394 -0.127488873976378;0.292454911850394 -0.131438433976378;0.291132561850394 -0.135351688976378;0.289771521850394 -0.139221218976378;0.288363576850394 -0.143048878976378;0.286911906850394 -0.146838113976378;0.285418631850394 -0.150583093976378;0.283886666850394 -0.154292033976378;0.282310711850394 -0.157956718976378;0.280693151850394 -0.161577413976378;0.279034251850394 -0.165164718976378;0.277334011850394 -0.168707768976378;0.275595081850394 -0.172212393976378;0.273817196850394 -0.175675148976378;0.271997971850394 -0.179096563976378;0.270139791850394 -0.182479288976378;0.268245571850394 -0.185820938976378;0.266309746850394 -0.189123368976378;0.264335231850394 -0.192384193976378;0.262324411850394 -0.195606328976378;0.260274636850394 -0.198787388976378;0.258186436850394 -0.201926843976378;0.256061666850394 -0.205027343976378;0.253903241850394 -0.208086503976378;0.251706126850394 -0.211109623976378;0.249472971850394 -0.214088488976378;0.247205896850394 -0.217028663976378;0.244902781850394 -0.219927498976378;0.242563096850394 -0.222790028976378;0.240192671850394 -0.225610953976378;0.237785941850394 -0.228393188976378;0.235342906850394 -0.231131433976378;0.232868866850394 -0.233833638976378;0.230361436850394 -0.236494238976378;0.227820086850394 -0.239115883976378;0.225245081850394 -0.241698838976378;0.222641721850394 -0.244243103976378;0.220002056850394 -0.246745763976378;0.217334036850394 -0.249209733976378;0.214635011850394 -0.251635013976378;0.211902066850394 -0.254021338976378;0.209138381850394 -0.256366058976378;0.206348461850394 -0.258674473976378;0.203525151850394 -0.260941813976378;0.200673221850394 -0.263167548976378;0.197792936850394 -0.265357243976378;0.194881646850394 -0.267505333976378;0.191944386850394 -0.269614733976378;0.188976121850394 -0.271687828976378;0.185981886850394 -0.273719583976378;0.182959296850394 -0.275712383976378;0.179908086850394 -0.277666493976378;0.176831171850394 -0.279581648976378;0.173730671850394 -0.281455463976378;0.170599431850394 -0.283295888976378;0.167442221850394 -0.285094708976378;0.164264341850394 -0.286854308976378;0.161055191850394 -0.288575483976378;0.157825636850394 -0.290260618976378;0.154570111850394 -0.291903883976378;0.151288616850394 -0.293508723976378;0.147986451850394 -0.295074343976378;0.144658581850394 -0.296601538976378;0.141307391850394 -0.298092163976378;0.137935266850394 -0.299544098976378;0.134537436850394 -0.300957343976378;0.131118936850394 -0.302331633976378;0.127676851850394 -0.303667233976378;0.124214361850394 -0.304963613976378;0.120730936850394 -0.306221568976378;0.117227106850394 -0.307443218976378;0.113702341850394 -0.308626178976378;0.110156906850394 -0.309772568976378;0.106593451850394 -0.310880268976378;0.103009061850394 -0.311946628976378;0.0994042668503938 -0.312976683976378;0.0957811868503937 -0.313970698976378;0.0921400868503938 -0.314922843976378;0.0884807018503938 -0.315841598976378;0.0848032968503938 -0.316721928976378;0.0811078718503938 -0.317560388976378;0.0773941618503938 -0.318365193976378;0.0736677318503938 -0.319128923976378;0.0699206318503938 -0.319861118976378;0.0661578968503938 -0.320549588976378;0.0623795268503937 -0.321199103976378;0.0585857868503937 -0.321815228976378;0.0547764118503937 -0.322389748976378;0.0509540518503938 -0.322930613976378;0.0471160568503938 -0.323432788976378;0.0432653418503938 -0.323898658976378;0.0393989918503938 -0.324325573976378;0.0355172718503938 -0.324714063976378;0.0316249518503938 -0.325065983976378;0.0277172618503938 -0.325381598976378;0.0237992368503938 -0.325661173976378;0.0198708768503937 -0.325899408976378;0.0159268818503938 -0.326103723976378;0.0119725518503937 -0.326269348976378;0.00800550185039377 -0.326401318976378;0.00402785185039376 -0.326494598976378;4.5166850393727e-05 -0.326548923976378;0.00070501685039379 -0.326548923976378;-0.00328270314960626 -0.326494598976378;-0.00725902814960622 -0.326401318976378;-0.0112250181496062 -0.326269348976378;-0.0151806731496063 -0.326103723976378;-0.0191233431496062 -0.325899408976378;-0.0230543531496062 -0.325661173976378;-0.0269723781496062 -0.325381598976378;-0.0308774181496062 -0.325065983976378;-0.0347707981496062 -0.324714063976378;-0.0386514581496063 -0.324325573976378;-0.0425178081496063 -0.323898658976378;-0.0463698481496062 -0.323432788976378;-0.0502075781496062 -0.322930613976378;-0.0540312631496062 -0.322389748976378;-0.0578406381496062 -0.321815228976378;-0.0616343781496062 -0.321199103976378;-0.0654127481496063 -0.320549588976378;-0.0691744231496063 -0.319861118976378;-0.0729199331496063 -0.319128923976378;-0.0766492781496062 -0.318365193976378;-0.0803613981496062 -0.317560388976378;-0.0840570881496062 -0.316721928976378;-0.0877344931496062 -0.315841598976378;-0.0913936131496062 -0.314922843976378;-0.0950349781496062 -0.313970698976378;-0.0986577931496063 -0.312976683976378;-0.102261528149606 -0.311946628976378;-0.105846978149606 -0.310880268976378;-0.109411758149606 -0.309772568976378;-0.112957193149606 -0.308626178976378;-0.116481958149606 -0.307443218976378;-0.119986053149606 -0.306221568976378;-0.123469213149606 -0.304963613976378;-0.126931968149606 -0.303667233976378;-0.130372463149606 -0.302331633976378;-0.133791228149606 -0.300957343976378;-0.137187733149606 -0.299544098976378;-0.140560918149606 -0.298092163976378;-0.143912373149606 -0.296601538976378;-0.147239978149606 -0.295074343976378;-0.150543468149606 -0.293508723976378;-0.153823638149606 -0.291903883976378;-0.157079163149606 -0.290260618976378;-0.160310308149606 -0.288575483976378;-0.163516543149606 -0.286854308976378;-0.166697073149606 -0.285094708976378;-0.169852958149606 -0.283295888976378;-0.172983138149606 -0.281455463976378;-0.176086023149606 -0.279581648976378;-0.179161613149606 -0.277666493976378;-0.182211498149606 -0.275712383976378;-0.185234088149606 -0.273719583976378;-0.188229648149606 -0.271687828976378;-0.191196588149606 -0.269614733976378;-0.194135173149606 -0.267505333976378;-0.197045403149606 -0.265357243976378;-0.199927013149606 -0.263167548976378;-0.202778943149606 -0.260941813976378;-0.205600928149606 -0.258674473976378;-0.208393233149606 -0.256366058976378;-0.211155858149606 -0.254021338976378;-0.213887213149606 -0.251635013976378;-0.216587828149606 -0.249209733976378;-0.219256908149606 -0.246745763976378;-0.221895248149606 -0.244243103976378;-0.224499933149606 -0.241698838976378;-0.227073613149606 -0.239115883976378;-0.229614963149606 -0.236494238976378;-0.232122658149606 -0.233833638976378;-0.234598023149606 -0.231131433976378;-0.237038408149606 -0.228393188976378;-0.239446463149606 -0.225610953976378;-0.241818213149606 -0.222790028976378;-0.244156308149606 -0.219927498976378;-0.246459688149606 -0.217028663976378;-0.248727823149606 -0.214088488976378;-0.250959918149606 -0.211109623976378;-0.253157033149606 -0.208086503976378;-0.255316518149606 -0.205027343976378;-0.257439963149606 -0.201926843976378;-0.259527103149606 -0.198787388976378;-0.261576613149606 -0.195606328976378;-0.263588758149606 -0.192384193976378;-0.265561948149606 -0.189123368976378;-0.267497773149606 -0.185820938976378;-0.269394643149606 -0.182479288976378;-0.271252823149606 -0.179096563976378;-0.273070723149606 -0.175675148976378;-0.274849933149606 -0.172212393976378;-0.276589128149606 -0.168707768976378;-0.278288043149606 -0.165164718976378;-0.279946943149606 -0.161577413976378;-0.281564238149606 -0.157956718976378;-0.283140193149606 -0.154292033976378;-0.284673748149606 -0.150583093976378;-0.286165698149606 -0.146838113976378;-0.287616043149606 -0.143048878976378;-0.289023988149606 -0.139221218976378;-0.290387678149606 -0.135351688976378;-0.291708438149606 -0.131438433976378;-0.292987063149606 -0.127488873976378;-0.294221433149606 -0.123495588976378;-0.295410488149606 -0.119465468976378;-0.296555553149606 -0.115389503976378;-0.297655568149606 -0.111276703976378;-0.298708678149606 -0.107120178976378;-0.299718063149606 -0.102922313976378;-0.300680808149606 -0.0986804589763782;-0.301598238149606 -0.0944020339763783;-0.302467703149606 -0.0900822689763783;-0.303290793149606 -0.0857190439763783;-0.304065653149606 -0.0813136839763783;-0.304794138149606 -0.0768698989763782;-0.305474923149606 -0.0723845089763782;-0.306106418149606 -0.0678580439763783;-0.306689948149606 -0.0632873239763783;-0.307224188149606 -0.0586776489763783;-0.307708078149606 -0.0540215989763783;-0.308144268149606 -0.0493295089763783;-0.308529843149606 -0.0445931639763783;-0.308865068149606 -0.0398178639763782;-0.309148353149606 -0.0349961889763783;-0.309382613149606 -0.0301355589763783;-0.309564933149606 -0.0252338539763782;-0.309695578149606 -0.0202905439763783;-0.309774548149606 -0.0153029789763782;-0.309800518149606 -0.0102743389763783;-0.309774548149606 -0.00521946397637828;-0.309695578149606 -0.000153193976378283;-0.309564933149606 0.00492341102362175;-0.309382613149606 0.0100106160236218;-0.309148353149606 0.0151052410236217;-0.308865068149606 0.0202104660236217;-0.308529843149606 0.0253194010236217;-0.308144268149606 0.0304376110236217;-0.307708078149606 0.0355608560236217;-0.307224188149606 0.0406891360236217;-0.306689948149606 0.0458216560236218;-0.306106418149606 0.0509565610236217;-0.305474923149606 0.0560935860236217;-0.304794138149606 0.0612311410236217;-0.304065653149606 0.0663713460236217;-0.303290793149606 0.0715112860236218;-0.302467703149606 0.0766501660236217;-0.301598238149606 0.0817877210236217;-0.300680808149606 0.0869226260236217;-0.299718063149606 0.0920548810236217;-0.298708678149606 0.0971821010236217;-0.297655568149606 0.102306671023622;-0.296555553149606 0.107424616023622;-0.295410488149606 0.112534876023622;-0.294221433149606 0.117638776023622;-0.292987063149606 0.122733666023622;-0.291708438149606 0.127821931023622;-0.290387678149606 0.132898536023622;-0.289023988149606 0.137964806023622;-0.287616043149606 0.143019681023622;-0.286165698149606 0.148062631023622;-0.284673748149606 0.153093656023622;-0.283140193149606 0.158109576023622;-0.281564238149606 0.163112511023622;-0.279946943149606 0.168099811023622;-0.278288043149606 0.173070416023622;-0.276589128149606 0.178024326023622;-0.274849933149606 0.182959951023622;-0.273070723149606 0.187878616023622;-0.271252823149606 0.192776876023622;-0.269394643149606 0.197655261023622;-0.267497773149606 0.202513241023622;-0.265561948149606 0.207347901023622;-0.263588758149606 0.212161891023622;-0.261576613149606 0.216951501023622;-0.259527103149606 0.221714876023622;-0.257439963149606 0.226456256023622;-0.255316518149606 0.231171931023622;-0.253157033149606 0.235860311023622;-0.250959918149606 0.240521661023622;-0.248727823149606 0.245153066023622;-0.246459688149606 0.249757441023622;-0.244156308149606 0.254332136023622;-0.241818213149606 0.258875561023622;-0.239446463149606 0.263387716023622;-0.237038408149606 0.267866746023622;-0.234598023149606 0.272314506023622;-0.232122658149606 0.276727551023622;-0.229614963149606 0.281105086023622;-0.227073613149606 0.285448171023622;-0.224499933149606 0.289754686023622;-0.221895248149606 0.294023836023622;-0.219256908149606 0.298255621023622;-0.216587828149606 0.302449511023622;-0.213887213149606 0.306603651023622;-0.211155858149606 0.310717246023622;-0.208393233149606 0.314788706023622;-0.205600928149606 0.318819886023622;-0.202778943149606 0.322808136023622;-0.199927013149606 0.326752396023622;-0.197045403149606 0.330652666023622;-0.194135173149606 0.334507886023622;-0.191196588149606 0.338316201023622;-0.188229648149606 0.342079201023622;-0.185234088149606 0.345794766023622;-0.182211498149606 0.349460776023622;-0.179161613149606 0.353079086023622;-0.176086023149606 0.356646781023622;-0.172983138149606 0.360162801023622;-0.169852958149606 0.363628471023622;-0.166697073149606 0.367040876023622;-0.163516543149606 0.370402931023622;-0.160310308149606 0.373709071023622;-0.157079163149606 0.376960091023622;-0.153823638149606 0.380156256023622;-0.150543468149606 0.383294386023622;-0.147239978149606 0.386376866023622;-0.143912373149606 0.389401046023622;-0.140560918149606 0.392367456023622;-0.137187733149606 0.395272386023622;-0.133791228149606 0.398119546023622;-0.130372463149606 0.400903106023622;-0.126931968149606 0.403624391023622;-0.123469213149606 0.406284991023622;-0.119986053149606 0.408880931023622;-0.116481958149606 0.411412211023622;-0.112957193149606 0.413878831023622;-0.109411758149606 0.416279201023622;-0.105846978149606 0.418612526023622;-0.102261528149606 0.420878276023622;-0.0986577931496063 0.423075656023622;-0.0950349781496062 0.425204401023622;-0.0913936131496062 0.427262126023622;-0.0877344931496062 0.429249891023622;-0.0840570881496062 0.431166106023622;-0.0803613981496062 0.433007591023622;-0.0766492781496062 0.434778056023622;-0.0729199331496063 0.436472996023622;-0.0691744231496063 0.438094796023622;-0.0654127481496063 0.439640011023622;-0.0616343781496062 0.441108641023622;-0.0578406381496062 0.442499891023622;-0.0540312631496062 0.443813231023622;-0.0502075781496062 0.445047866023622;-0.0463698481496062 0.446203266023622;-0.0425178081496063 0.447277576023622;-0.0386514581496063 0.448270001023622;-0.0347707981496062 0.449181071023622;-0.0308774181496062 0.450008136023622;-0.0269723781496062 0.450753316023622;-0.0230543531496062 0.451413166023622;-0.0191233431496062 0.451986626023622;-0.0151806731496063 0.452475551023622;-0.0112250181496062 0.452875701023622;-0.00725902814960622 0.453190256023622;-0.00328270314960626 0.453413916023622;0.00070501685039379 0.453549861023622];
outline{3} = 5.5*[-0.0518871481496062 0.444574046023622;-0.0503384881496063 0.446502186023622;-0.0487895631496062 0.448430591023622;-0.0472019481496062 0.450366416023622;-0.0456220181496062 0.452268586023622;-0.0440291031496063 0.454155386023622;-0.0424258531496062 0.456026816023622;-0.0407889481496063 0.457921301023622;-0.0391650281496062 0.459779481023622;-0.0375257381496062 0.461622291023622;-0.0358758481496062 0.463444431023622;-0.0342145631496062 0.465261006023622;-0.0325453281496062 0.467072811023622;-0.0308360781496062 0.468902901023622;-0.0291125181496062 0.470714441023622;-0.0274186381496063 0.472492591023622;-0.0256937531496062 0.474244771023622;-0.0239715181496062 0.475999336023622;-0.0222325881496062 0.477736146023622;-0.0204494031496062 0.479490976023622;-0.0186691331496062 0.481237856023622;-0.0168846231496063 0.482940746023622;-0.0151067381496062 0.484633566023622;-0.0133081831496062 0.486308101023622;-0.0107398031496062 0.488608831023622;-0.00906023314960622 0.489802126023622;-0.00681647814960622 0.491008141023622;-0.00452104814960624 0.491755971023622;-0.00226192314960626 0.492356726023622;1.41618503937522e-05 0.492656971023622;0.000733371850393752 0.492656971023622;0.0030081318503938 0.492356726023622;0.00526725685039378 0.491755971023622;0.00756295185039377 0.491008141023622;0.00980644185039375 0.489802126023622;0.0114860118503938 0.488608831023622;0.0140559818503938 0.486308101023622;0.0158518868503938 0.484633566023622;0.0176324218503938 0.482940746023622;0.0194153418503938 0.481237856023622;0.0211958768503937 0.479490976023622;0.0229787968503937 0.477736146023622;0.0247179918503938 0.475999336023622;0.0264415518503938 0.474244771023622;0.0281651118503937 0.472492591023622;0.0298574018503938 0.470714441023622;0.0315809618503938 0.468902901023622;0.0332915368503938 0.467072811023622;0.0349607718503937 0.465261006023622;0.0366223218503938 0.463444431023622;0.0382732718503937 0.461622291023622;0.0399115018503938 0.459779481023622;0.0415340968503938 0.457921301023622;0.0431723268503938 0.456026816023622;0.0447742518503938 0.454155386023622;0.0463682268503937 0.452268586023622;0.0479468318503938 0.450366416023622;0.0495357718503937 0.448430591023622;0.0510833718503938 0.446502186023622;0.0526336218503938 0.444574046023622];
outline{4} = 5.5*[-0.304627188149606 0.0624077410236217;-0.305125388149606 0.0617839310236217;-0.306058453149606 0.0605803010236218;-0.307023583149606 0.0593032660236218;-0.308017333149606 0.0579546810236217;-0.309038378149606 0.0565390510236218;-0.310080093149606 0.0550598210236218;-0.311142213149606 0.0535196410236218;-0.312220233149606 0.0519216910236217;-0.313312298149606 0.0502691510236217;-0.314414698149606 0.0485649360236218;-0.315525048149606 0.0468127560236218;-0.316640433149606 0.0450165860236217;-0.317755553149606 0.0431777510236217;-0.318870938149606 0.0412999610236218;-0.319981288149606 0.0393885160236217;-0.321083688149606 0.0374447410236217;-0.322177078149606 0.0354715510236217;-0.323254833149606 0.0334721260236217;-0.324318543149606 0.0314507060236218;-0.325361583149606 0.0294099410236217;-0.326382363149606 0.0273538060236217;-0.327378763149606 0.0252844210236218;-0.328345218149606 0.0232062910236217;-0.329280933149606 0.0211225960236217;-0.330181403149606 0.0190338660236217;-0.331045833149606 0.0169456660236217;-0.331868658149606 0.0148619710236217;-0.332649083149606 0.0127835760236217;-0.333381278149606 0.0107157810236218;-0.334065773149606 0.0086607060236217;-0.334697268149606 0.00662126602362171;-0.335273113149606 0.00460117102362175;-0.335791983149606 0.00260439602362175;-0.336248843149606 0.000633591023621698;-0.336640778149606 -0.00130859397637823;-0.336965668149606 -0.0032187139763783;-0.337219273149606 -0.00509411897637826;-0.337400268149606 -0.00693030397637829;-0.337505208149606 -0.00872382397637825;-0.337529853149606 -0.0104733539763783;-0.337472878149606 -0.0121741239763782;-0.337330573149606 -0.0138224239763783;-0.337100288149606 -0.0154182539763783;-0.336778048149606 -0.0169541939763783;-0.336360143149606 -0.0184294489763782;-0.335847633149606 -0.0198397789763783;-0.335233098149606 -0.0211833289763783;-0.334514948149606 -0.0224563889763783;-0.333690533149606 -0.0236549839763783;-0.332757733149606 -0.0247754039763782;-0.331710983149606 -0.0258160589763783;-0.330548958149606 -0.0267735039763783;-0.329269273149606 -0.0276432339763782;-0.327869278149606 -0.0284249839763783;-0.326343673149606 -0.0291108039763783;-0.324689808149606 -0.0297006939763783;-0.322908213149606 -0.0301925339763783;-0.320991733149606 -0.0305810239763783;-0.318938248149606 -0.0308629839763782;-0.316746433149606 -0.0310339089763783;-0.314412048149606 -0.0310959189763783;-0.311931648149606 -0.0310389439763783;-0.309303643149606 -0.0308629839763782]; 
outline{5} = 5.5*[0.305376311850394 0.0624064160236217;0.305873186850394 0.0617839310236217;0.306804661850394 0.0605803010236218;0.307770056850394 0.0593032660236218;0.308763806850394 0.0579546810236217;0.309785911850394 0.0565390510236218;0.310826301850394 0.0550598210236218;0.311887361850394 0.0535196410236218;0.312966441850394 0.0519216910236217;0.314058506850394 0.0502691510236217;0.315161171850394 0.0485649360236218;0.316271256850394 0.0468127560236218;0.317386641850394 0.0450165860236217;0.318502026850394 0.0431777510236217;0.319617411850394 0.0412999610236218;0.320727496850394 0.0393885160236217;0.321830161850394 0.0374447410236217;0.322922226850394 0.0354715510236217;0.324001306850394 0.0334721260236217;0.325065016850394 0.0314507060236218;0.326107791850394 0.0294099410236217;0.327130161850394 0.0273538060236217;0.328123911850394 0.0252844210236218;0.329091691850394 0.0232062910236217;0.330028466850394 0.0211225960236217;0.330928936850394 0.0190338660236217;0.331793366850394 0.0169456660236217;0.332616456850394 0.0148619710236217;0.333395291850394 0.0127835760236217;0.334127751850394 0.0107157810236218;0.334810921850394 0.0086607060236217;0.335442416850394 0.00662126602362171;0.336019321850394 0.00460117102362175;0.336539516850394 0.00260439602362175;0.336995051850394 0.000633591023621698;0.337388311850394 -0.00130859397637823;0.337711876850394 -0.0032187139763783;0.337965481850394 -0.00509411897637826;0.338146741850394 -0.00693030397637829;0.338252741850394 -0.00872382397637825;0.338276061850394 -0.0104733539763783;0.338219086850394 -0.0121741239763782;0.338076781850394 -0.0138224239763783;0.337846496850394 -0.0154182539763783;0.337522931850394 -0.0169541939763783;0.337106351850394 -0.0184294489763782;0.336593841850394 -0.0198397789763783;0.335980631850394 -0.0211833289763783;0.335261156850394 -0.0224563889763783;0.334438331850394 -0.0236549839763783;0.333503941850394 -0.0247754039763782;0.332455866850394 -0.0258160589763783;0.331296491850394 -0.0267735039763783;0.330015481850394 -0.0276432339763782;0.328615486850394 -0.0284249839763783;0.327088556850394 -0.0291108039763783;0.325437606850394 -0.0297006939763783;0.323654421850394 -0.0301925339763783;0.321736881850394 -0.0305810239763783;0.319684721850394 -0.0308629839763782;0.317492641850394 -0.0310339089763783;0.315158521850394 -0.0310959189763783;0.312679181850394 -0.0310389439763783;0.310049851850394 -0.0308629839763782];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate a square outline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_square(layout)

% get index of all relevant channels
ind1 = strcmp(layout.label,'COMNT');
ind2 = strcmp(layout.label,'SCALE');
ind  = ~(ind1 | ind2);

x = layout.pos(ind,1);
y = layout.pos(ind,2);
w = layout.width(ind);
h = layout.height(ind);

% determine the bounding box and add a little space around it
space = min(mean(w), mean(h))/4;
xmin = min(x - 0.5*w - space);
xmax = max(x + 0.5*w + space);
ymin = min(y - 0.5*h - space);
ymax = max(y + 0.5*h + space);

% construct the outline
outline = {[
  xmin ymax
  xmax ymax
  xmax ymin
  xmin ymin
  xmin ymax
  ]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate an outline from the boundary/convex hull of pos+width/height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_convex(layout)

% get index of all relevant channels
ind1 = strcmp(layout.label,'COMNT');
ind2 = strcmp(layout.label,'SCALE');
ind  = ~(ind1 | ind2);
% get position of all channel 'boxes' and draw boundary around it, or, in older matlabs, a convex hull
x = layout.pos(ind,1);
y = layout.pos(ind,2);
w = layout.width(ind);
h = layout.height(ind);
boxpos = [
  x - (w/2) y - (h/2);  % lb
  x - (w/2) y + (h/2);  % lt
  x + (w/2) y - (h/2);  % rb
  x + (w/2) y + (h/2);  % rt
  ];
if ft_platform_supports('boundary')
  k = boundary(boxpos,.2);
  outline{1} = boxpos(k,:);
else
  outline{1} = boxpos(convhull(boxpos(:,1),boxpos(:,2)),:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate an outline from the headshape or anatomical mri
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_headshape(cfg, sens)

if ~isempty(cfg.headshape)
  if ischar(cfg.headshape) && exist(cfg.headshape, 'file')
    ft_info('reading headshape from file %s\n', cfg.headshape);
    outlbase = ft_read_headshape(cfg.headshape);
  elseif isstruct(cfg.headshape)
    outlbase = cfg.headshape;
  else
    ft_error('incorrect specification of cfg.headshape')
  end
elseif ~isempty(cfg.mri)
  if ischar(cfg.mri) && exist(cfg.mri, 'file')
    ft_info('reading MRI from file %s\n', cfg.mri);
    outlbase = ft_read_mri(cfg.mri);
  elseif ft_datatype(cfg.mri, 'volume')
    outlbase = cfg.mri;
  else
    ft_error('incorrect specification of cfg.mri')
  end
  % create mesh from anatomical field, this will be used as headshape below
  cfgpm = [];
  cfgpm.method      = 'projectmesh';
  cfgpm.tissue      = 'brain';
  cfgpm.numvertices = 1e5;
  outlbase = ft_prepare_mesh(cfgpm, outlbase);
end

% check that we have the right data in outlbase
assert(isfield(outlbase, 'pos'), 'the headshape does not contain any vertices')

% check coordinate system of outlbase
assert(isfield(outlbase, 'coordsys'), 'no coordsys field found in headshape/mri, use ft_determine_coordsys')
assert(isfield(sens, 'coordsys'), 'no coordsys field found in sensor structure, use ft_determine_coordsys')
assert(isequal(outlbase.coordsys, sens.coordsys), 'the coordinate system of headshape/mri does not match that of sensors')

% match head geometry units with that of the sensors
outlbase = ft_convert_units(outlbase, sens.unit);

% there can be multiple meshes, e.g. left and right hemispheres
outline = cell(size(outlbase));
for i=1:numel(outlbase)
  % generate outline based on matlab version
  if ft_platform_supports('boundary')
    
    % extract points indicating brain
    braincoords = outlbase(i).pos;
    % apply projection and extract XY
    if isequal(cfg.projection,'orthographic') && ~isempty(cfg.viewpoint)
      braincoords = getorthoviewpos(braincoords, outlbase(i).coordsys, cfg.viewpoint);
    else
      % project identically as in sens2lay using cfg.rotate and elproj (this should be kept identical to sens2lay)
      if isempty(cfg.rotate)
        switch ft_senstype(sens)
          case {'ctf151', 'ctf275', 'bti148', 'bti248', 'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar', 'yokogawa160', 'yokogawa160_planar', 'yokogawa64', 'yokogawa64_planar', 'yokogawa440', 'yokogawa440_planar', 'magnetometer', 'meg'}
            rotatez = 90;
          case {'neuromag122', 'neuromag306'}
            rotatez = 0;
          case 'electrode'
            rotatez = 90;
          otherwise
            rotatez = 0;
        end
      end
      braincoords = ft_warp_apply(rotate([0 0 rotatez]), braincoords, 'homogenous');
      braincoords = elproj(braincoords, cfg.projection);
    end
    
    % get outline
    k = boundary(braincoords,.8);
    outline{i} = braincoords(k,:);
    
  else % matlab version fallback
    
    % plot mesh in rotated view, rotate, screencap, and trace frame to generate outline
    if isequal(cfg.projection,'orthographic') && ~isempty(cfg.viewpoint)
      outlbase(i).pos = getorthoviewpos(outlbase(i).pos, outlbase(i).coordsys, cfg.viewpoint);
    else
      % project identically as in sens2lay using cfg.rotate and elproj (this should be kept identical to sens2lay)
      if isempty(cfg.rotate)
        switch ft_senstype(sens)
          case {'ctf151', 'ctf275', 'bti148', 'bti248', 'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar', 'yokogawa160', 'yokogawa160_planar', 'yokogawa64', 'yokogawa64_planar', 'yokogawa440', 'yokogawa440_planar', 'magnetometer', 'meg'}
            rotatez = 90;
          case {'neuromag122', 'neuromag306'}
            rotatez = 0;
          case 'electrode'
            rotatez = 90;
          otherwise
            rotatez = 0;
        end
      end
      outlbase(i).pos = ft_warp_apply(rotate([0 0 rotatez]), outlbase(i).pos, 'homogenous');
      outlbase(i).pos = elproj(outlbase(i).pos, cfg.projection);
    end
    h = figure('visible', 'off');
    ft_plot_mesh(outlbase(i), 'facecolor', [0 0 0], 'EdgeColor', 'none');
    view([0 90])
    axis tight
    xlim = get(gca, 'xlim');
    ylim = get(gca, 'ylim');
    set(gca,'OuterPosition', [0.2 0.2 .6 .6]) % circumvent weird matlab bug, wth?
    % extract frame for tracing
    drawnow % need to flush buffer, otherwise frame will not get extracted properly
    frame     = getframe(h);
    close(h)
    imtotrace = double(~logical(sum(frame.cdata,3))); % needs to be binary to trace
    % image is not in regular xy space, flip, and transpose
    imtotrace = flipud(imtotrace).';
    
    % trace image generated above
    [row, col] = find(imtotrace, 1, 'first'); % set an arbitrary starting point
    trace = bwtraceboundary(imtotrace, [row, col], 'N');
    
    % convert to sens coordinates
    x = trace(:,1);
    y = trace(:,2);
    x = x - min(x);
    x = x ./ max(x);
    x = x .* (xlim(2)-xlim(1));
    x = x - abs(xlim(1));
    y = y - min(y);
    y = y ./ max(y);
    y = y .* (ylim(2)-ylim(1));
    y = y - abs(ylim(1));
    
    outline{i} = [x y];
  end
  
  % subsample the outline
  if size(outline{i},1)>5e3 % 5e3 points should be more than enough to get an outine with acceptable detail
    outline{i} = outline{i}(1:floor(size(outline,1)/5e3):end,:);
  end
  
end % for numel(outlbase)
