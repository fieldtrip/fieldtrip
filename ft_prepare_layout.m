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
%   cfg.outline     = string, how to create the outline, can be 'circle', 'square', 'convex', 'headshape', 'mri' or 'no' (default is automatic)
%   cfg.mask        = string, how to create the mask, can be 'circle', 'square', 'convex', 'headshape', 'mri' or 'no' (default is automatic)
%   cfg.headshape   = surface mesh (e.g. pial, head, etc) to be used for generating an outline, see FT_READ_HEADSHAPE for details
%   cfg.mri         = segmented anatomical MRI to be used for generating an outline, see FT_READ_MRI and FT_VOLUMESEGMENT for details
%   cfg.montage     = 'no' or a montage structure (default = 'no')
%   cfg.image       = filename, use an image to construct a layout (e.g. useful for ECoG grids)
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

if isa(cfg.layout, 'config')
  % convert the nested config-object back into a normal structure
  cfg.layout = struct(cfg.layout);
end

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
  layout.pos     = zeros(nchan,2);  % centered at (0,0)
  layout.width   = ones(nchan,1) * 1.0;
  layout.height  = ones(nchan,1) * 1.0;
  layout.mask    = {};
  layout.outline = {};
  
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
  
  if isempty(strfind(cfg.layout, '.'))
    
    % check whether a corresponding mat or lay exists
    if exist([cfg.layout '.mat'], 'file')
      ft_info('appending .mat to layout file\n');
      cfg.layout = [cfg.layout '.mat'];
      layout = ft_prepare_layout(cfg);
      return
    else
      ft_info('appending .lay to layout file\n');
      cfg.layout = [cfg.layout '.lay'];
      layout = ft_prepare_layout(cfg);
      return
    end
    
  elseif ft_filetype(cfg.layout, 'matlab')
    
    ft_info('reading layout from file %s\n', cfg.layout);
    if ~exist(cfg.layout, 'file')
      ft_error('the specified layout file %s was not found', cfg.layout);
    end
    tmp = load(cfg.layout, 'lay*');
    if isfield(tmp, 'layout')
      layout = tmp.layout;
    elseif isfield(tmp, 'lay')
      layout = tmp.lay;
    else
      ft_error('mat file does not contain a layout');
    end
    
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
    
    % assume that cfg.layout is an electrode file
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
    ft_info('reading background image from %s\n', cfg.image);
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
  
  if strcmp(cfg.outline, 'circle') || strcmp(cfg.mask, 'circle')
    % Scale the electrode positions to fit within a unit circle, i.e. electrode radius = 0.45
    ind_scale = find(strcmp('SCALE', layout.label));
    ind_comnt = find(strcmp('COMNT', layout.label));
    sel = setdiff(1:length(layout.label), [ind_scale ind_comnt]); % these are excluded for scaling
    x = layout.pos(sel,1);
    y = layout.pos(sel,2);
    if istrue(cfg.center)
      % the following centers all electrodes around zero
      xrange = range(x);
      yrange = range(y);
      shiftx = min(x);
      shifty = min(y);
    else
      % the following prevent topography distortion in case electrodes are not evenly distributed over the whole head
      xrange = 2*( max(max(x),abs(min(x)) ));
      yrange = 2*( max(max(y),abs(min(y)) ));
      shiftx = ( max(max(x),abs(min(x)) )).*sign(min(x));
      shifty = ( max(max(y),abs(min(y)) )).*sign(min(y));
    end
    if xrange==0
      xrange = 1;
    end
    if yrange==0
      yrange = 1;
    end
    % First scale the width and height of the box for multiplotting
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
      case 'convex'
        layout.outline = outline_convex(layout);
      case 'square'
        layout.outline = outline_square(layout);
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
        layout.mask = layout.mask(1); % the first is the circle, the others are nose and ears
      case 'convex'
        layout.mask = outline_convex(layout);
      case 'square'
        layout.mask = outline_square(layout);
      case {'headshape', 'mri'}
        % the configuration should contain the headshape or mri
        % the (segmented) mri will be converted into a headshape on the fly
        if isequal(cfg.mask,cfg.outline) && exist('hsoutline','var')
          layout.mask = hsoutline;
        else
          layout.mask = outline_headshape(cfg, sens);
        end
      otherwise
        layout.mask = {};
    end
  end
  
end % create outline if style=2d


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
