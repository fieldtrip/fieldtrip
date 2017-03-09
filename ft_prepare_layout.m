function [layout, cfg] = ft_prepare_layout(cfg, data)

% FT_PREPARE_LAYOUT loads or creates a 2-D layout of the channel locations.
% This layout is required for plotting the topographical distribution of
% the potential or field distribution, or for plotting timecourses in a
% topographical arrangement.
%
% Use as
%   layout = ft_prepare_layout(cfg, data)
%
% There are several ways in which a 2-D layout can be made: it can be read
% directly from a *.mat file containing a variable 'lay', it can be created
% based on 3-D electrode or gradiometer positions in the configuration or
% in the data, or it can be created based on the specification of an
% electrode or gradiometer file. Layouts can also come from an ASCII *.lay
% file, but this type of layout is no longer recommended.
%
% You can specify any one of the following configuration options
%   cfg.layout      = filename containg the layout (.mat or .lay file)
%                     can also be a layout structure, which is simply
%                     returned as-is (see below for details)
%   cfg.rotate      = number, rotation around the z-axis in degrees (default = [], which means automatic)
%   cfg.projection  = string, 2D projection method can be 'stereographic', 'orthographic',
%                     'polar', 'gnomic' or 'inverse' (default = 'polar')
%                     When 'orthographic', cfg.viewpoint can be used to indicate to specificy projection (keep empty for legacy projection)
%   cfg.viewpoint   = string indicating 'viewpoint' used for orthographic projection of 3D electrode coordinates to 2D plane
%                     (requires cfg.projection = 'orthographic')  
%                     Viewpoints are as follows:
%                     'left'      = left  sagittal view,    L=anterior, R=posterior, top=top, bottom=bottom
%                     'right'     = right sagittal view,    L=posterior, R=anterior, top=top, bottom=bottom
%                     'inferior'  = inferior axial view,    L=R, R=L, top=anterior, bottom=posterior
%                     'superior'  = superior axial view,    L=L, R=R, top=anterior, bottom=posterior
%                     'anterior'  = anterior  coronal view, L=R, R=L, top=top, bottom=bottom
%                     'posterior' = posterior coronal view, L=L, R=R, top=top, bottom=bottom
%                     'auto'      = automatic guess of the most optimal of the above
%                      tip: use cfg.viewpoint = auto per iEEG electrode grid/strip/depth for more accurate results
%                      tip: to obtain overview of e.g. all iEEG electrodes, choose superior/inferior, use cfg.headshape/mri, and
%                           plot using ft_layoutplot with cfg.box/mask = 'no'
%   cfg.headshape   = surface mesh (e.g. pial, head, etc) to be used for generating a layout outline using cfg.viewpoint
%                     If used, needs to be in same coordinate space as cfg.elec/grad/opto/file. See FT_READ_HEADSHAPE
%   cfg.mri         = mri to be used for generating a brain outline as layout outline for cfg.viewpoint
%                     If used, needs to be in same coordinate space as cfg.elec/grad/opto/file. 
%                     Mri needs to be segmented and contain a 'brain' field. If not, segmentation is attempted automatically.
%                     See FT_READ_MRI and FT_VOLUMESEGMENT.
%   cfg.elec        = structure with electrode definition, or
%   cfg.elecfile    = filename containing electrode definition
%   cfg.grad        = structure with gradiometer definition, or
%   cfg.gradfile    = filename containing gradiometer definition
%   cfg.opto        = structure with optode structure definition, or
%   cfg.optofile    = filename containing optode structure definition
%   cfg.output      = filename (ending in .mat or .lay) to which the layout
%                     will be written (default = [])
%   cfg.montage     = 'no' or a montage structure (default = 'no')
%   cfg.image       = filename, use an image to construct a layout (e.g. useful for ECoG grids)
%   cfg.bw          = if an image is used and bw = 1 transforms the image in
%                     black and white (default = 0, do not transform)
%   cfg.overlap     = string, how to deal with overlapping channels when
%                     layout is constructed from a sensor configuration
%                     structure (can be 'shift' (shift the positions in 2D
%                     space to remove the overlap (default)), 'keep' (don't
%                     shift, retain the overlap), 'no' (throw error when
%                     overlap is present))
%   cfg.channel     = 'all', or Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%                     (can be used with most ways for generating a layout)
%   cfg.boxchannel  = 'all', or Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%                      specificies channels to use for determining channel box size (default = all, recommended for MEG/EEG, selection recommended for iEEG)
%   cfg.skipscale   = 'yes' or 'no', whether the scale should be included in the layout or not (default = 'no')
%   cfg.skipcomnt   = 'yes' or 'no', whether the comment should be included in the layout or not (default = 'no')
%
% Alternatively the layout can be constructed from either
%   data.elec     = structure with electrode positions
%   data.grad     = structure with gradiometer definition
%   data.opto     = structure with optode structure definition
%
% Alternatively you can specify the following layouts which will be
% generated for all channels present in the data. Note that these layouts
% are suitable for multiplotting, but not for topoplotting.
%   cfg.layout = 'ordered'    will give you a NxN ordered layout
%   cfg.layout = 'vertical'   will give you a Nx1 ordered layout
%   cfg.layout = 'horizontal' will give you a 1xN ordered layout
%   cfg.layout = 'butterfly'  will give you a layout with all channels on top of each other
%   cfg.layout = 'circular'   will distribute the channels on a circle
%
% The output layout structure will contain the following fields
%   layout.label   = Nx1 cell-array with channel labels
%   layout.pos     = Nx2 matrix with channel positions
%   layout.width   = Nx1 vector with the width of each box for multiplotting
%   layout.height  = Nx1 matrix with the height of each box for multiplotting
%   layout.mask    = optional cell-array with line segments that determine the area for topographic interpolation
%   layout.outline = optional cell-array with line segments that represent
%                    the head, nose, ears, sulci or other anatomical features
%
% See also FT_TOPOPLOTER, FT_TOPOPLOTTFR, FT_MULTIPLOTER, FT_MULTIPLOTTFR, FT_PLOT_LAY

% undocumented and non-recommended option (for SPM only)
%   cfg.style       string, '2d' or '3d' (default = '2d')

% Copyright (C) 2007-2013, Robert Oostenveld
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
hasdata = exist('data', 'var');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic check/initialization of input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~hasdata
  data = struct([]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default configuration options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.rotate       = ft_getopt(cfg, 'rotate',     []); % [] => rotation is determined based on the type of sensors
cfg.style        = ft_getopt(cfg, 'style',      '2d');
cfg.projection   = ft_getopt(cfg, 'projection', 'polar');
cfg.layout       = ft_getopt(cfg, 'layout',     []);
cfg.grad         = ft_getopt(cfg, 'grad',       []);
cfg.elec         = ft_getopt(cfg, 'elec',       []);
cfg.opto         = ft_getopt(cfg, 'opto',       []);
cfg.gradfile     = ft_getopt(cfg, 'gradfile',   []);
cfg.elecfile     = ft_getopt(cfg, 'elecfile',   []);
cfg.optofile     = ft_getopt(cfg, 'optofile',   []);
cfg.output       = ft_getopt(cfg, 'output',     []);
cfg.feedback     = ft_getopt(cfg, 'feedback',   'no');
cfg.montage      = ft_getopt(cfg, 'montage',    'no');
cfg.image        = ft_getopt(cfg, 'image',      []);
cfg.mesh         = ft_getopt(cfg, 'mesh',       []); % experimental, should only work with meshes defined in 2D
cfg.bw           = ft_getopt(cfg, 'bw',         0);
cfg.channel      = ft_getopt(cfg, 'channel',    'all');
cfg.skipscale    = ft_getopt(cfg, 'skipscale',  'no');
cfg.skipcomnt    = ft_getopt(cfg, 'skipcomnt',  'no');
cfg.boxchannel   = ft_getopt(cfg, 'boxchannel', 'all');
cfg.overlap      = ft_getopt(cfg, 'overlap',    'shift');
cfg.viewpoint    = ft_getopt(cfg, 'viewpoint',  []);
cfg.headshape    = ft_getopt(cfg, 'headshape',  []); % separate form cfg.mesh
cfg.mri          = ft_getopt(cfg, 'mri',        []);

% headshape/mri are mutually exclusive
if ~isempty(cfg.headshape) && ~isempty(cfg.mri)
  error('cfg.headshape and cfg.mri are mutually exclusive, please use only one of the two')
end
% cfg.viewpoint can only be used together with cfg.projection = 'orthographic'
if ~isempty(cfg.viewpoint) && ~isequal(cfg.projection,'orthographic')
  error('cfg.viewpoint can only used in the case of orthographic projection')
end
if ~isempty(cfg.viewpoint) && ~isempty(cfg.rotate)
  error('cfg.viewpoint and cfg.rotate cannot be used jointly')
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
      error('the number of elements in the polar angle vector should be equal to the number of channels');
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
  skipscale = true; % a scale is not desired
  skipcomnt = true; % a comment is initially not desired, or at least requires more thought
  
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
  skipscale = true; % a scale is not desired
  skipcomnt = true; % a comment is initially not desired, or at least requires more thought
  
elseif isequal(cfg.layout, 'vertical') || isequal(cfg.layout,'horizontal')
  if hasdata && ~isempty(data)
    % look at the data to determine the overlapping channels
    data = ft_checkdata(data);
    originalorder = cfg.channel;
    cfg.channel = ft_channelselection(cfg.channel, data.label);
    if iscell(originalorder) && length(originalorder)==length(cfg.channel)
      % try to keep the order identical to that specified in the configuration
      [sel1, sel2] = match_str(originalorder, cfg.channel);
      % re-order them according to the cfg specified by the user
      cfg.channel  = cfg.channel(sel2);
    end
    assert(iscell(cfg.channel), 'cfg.channel should be a valid set of channels');
    nchan        = length(cfg.channel);
    layout.label = cfg.channel;
  else
    assert(iscell(cfg.channel), 'cfg.channel should be a valid set of channels');
    nchan        = length(cfg.channel);
    layout.label = cfg.channel;
  end
  for i=1:(nchan+2)
    switch cfg.layout
      case 'vertical'
        x = 0.5;
        y = 1-i/(nchan+1+2);
        layout.pos   (i,:) = [x y];
        layout.width (i,1) = 0.9;
        layout.height(i,1) = 0.9 * 1/(nchan+1+2);
      case 'horizontal'
        x = i/(nchan+1+2);
        y = 0.5;
        layout.pos   (i,:) = [x y];
        layout.width (i,1) = 0.9 * 1/(nchan+1+2);
        layout.height(i,1) = 0.9;
    end
    if i==(nchan+1)
      layout.label{i}   = 'SCALE';
    elseif i==(nchan+2)
      layout.label{i}   = 'COMNT';
    end
  end
  layout.mask    = {};
  layout.outline = {};
  
elseif any(strcmp(cfg.layout, {'1column', '2column', '3column', '4column', '5column', '6column', '7column', '8column', '9column'}))
  % it can be 2column, 3column, etcetera
  % note that this code (in combination with the code further down) fails for 1column
  if hasdata && ~isempty(data)
    % look at the data to determine the overlapping channels
    data = ft_checkdata(data);
    originalorder = cfg.channel;
    cfg.channel = ft_channelselection(cfg.channel, data.label);
    if iscell(originalorder) && length(originalorder)==length(cfg.channel)
      % try to keep the order identical to that specified in the configuration
      [sel1, sel2] = match_str(originalorder, cfg.channel);
      % re-order them according to the cfg specified by the user
      cfg.channel  = cfg.channel(sel2);
    end
    assert(iscell(cfg.channel), 'cfg.channel should be a valid set of channels');
    nchan        = length(cfg.channel);
    layout.label = cfg.channel;
  else
    assert(iscell(cfg.channel), 'cfg.channel should be a valid set of channels');
    nchan        = length(cfg.channel);
    layout.label = cfg.channel;
  end
  
  ncol = find(strcmp(cfg.layout, {'1column', '2column', '3column', '4column', '5column', '6column', '7column', '8column', '9column'}));
  nrow = ceil(nchan/ncol);
  
  k = 0;
  for i=1:ncol
    for j=1:nrow
      k = k+1;
      if k>nchan
        continue
      end
      x = i/ncol - 1/(ncol*2);
      y = 1-j/(nrow+1);
      layout.pos   (k,:) = [x y];
      layout.width (k,1) = 0.85/ncol;
      layout.height(k,1) = 0.9 * 1/(nrow+1);
    end
  end
  
  layout.mask    = {};
  layout.outline = {};
  skipscale = true; % a scale is not desired
  skipcomnt = true; % a comment is initially not desired, or at least requires more thought
  
elseif isequal(cfg.layout, 'ordered')
  if hasdata && ~isempty(data)
    % look at the data to determine the overlapping channels
    data = ft_checkdata(data);
    cfg.channel = ft_channelselection(cfg.channel, data.label);
    chanindx    = match_str(data.label, cfg.channel);
    nchan       = length(data.label(chanindx));
    layout.label   = data.label(chanindx);
  else
    assert(iscell(cfg.channel), 'cfg.channel should be a valid set of channels');
    nchan        = length(cfg.channel);
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
        layout.pos(k,:)    = [x y];
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
      return;
    else
      cfg.layout = [cfg.layout(1:end-3) 'lay'];
      layout = ft_prepare_layout(cfg);
      return;
    end
    
  elseif ft_filetype(cfg.layout, 'matlab')
    
    fprintf('reading layout from file %s\n', cfg.layout);
    if ~exist(cfg.layout, 'file')
      error('the specified layout file %s was not found', cfg.layout);
    end
    tmp = load(cfg.layout, 'lay*');
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
      ft_warning(sprintf('layout file %s was not found on your path, attempting to use a similarly named .mat file instead',cfg.layout));
      cfg.layout = [cfg.layout(1:end-3) 'mat'];
      layout = ft_prepare_layout(cfg);
      return;
    end
    
  elseif ~ft_filetype(cfg.layout, 'layout')
    % assume that cfg.layout is an electrode file
    fprintf('creating layout from electrode file %s\n', cfg.layout);
    sens = ft_read_sens(cfg.layout);
    layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  end
  
elseif ischar(cfg.elecfile)
  fprintf('creating layout from electrode file %s\n', cfg.elecfile);
  sens = ft_read_sens(cfg.elecfile);
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  
elseif ~isempty(cfg.elec) && isstruct(cfg.elec)
  fprintf('creating layout from cfg.elec\n');
  sens = ft_datatype_sens(cfg.elec);
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  
elseif isfield(data, 'elec') && isstruct(data.elec)
  fprintf('creating layout from data.elec\n');
  data = ft_checkdata(data);
  sens = ft_datatype_sens(data.elec);
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  
elseif ischar(cfg.gradfile)
  fprintf('creating layout from gradiometer file %s\n', cfg.gradfile);
  sens = ft_read_sens(cfg.gradfile);
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  
elseif ~isempty(cfg.grad) && isstruct(cfg.grad)
  fprintf('creating layout from cfg.grad\n');
  sens = ft_datatype_sens(cfg.grad);
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
  
elseif isfield(data, 'grad') && isstruct(data.grad)
  fprintf('creating layout from data.grad\n');
  data = ft_checkdata(data);
  sens = ft_datatype_sens(data.grad);
  layout = sens2lay(sens, cfg.rotate, cfg.projection, cfg.style, cfg.overlap, cfg.viewpoint, cfg.boxchannel);
 
elseif ischar(cfg.optofile)
  fprintf('creating layout from optode file %s\n', cfg.optofile);
  sens = ft_read_sens(cfg.optofile);
  if (hasdata)
    layout = opto2lay(sens, data.label);
  else
    layout = opto2lay(sens, sens.label);
  end

elseif ~isempty(cfg.opto) && isstruct(cfg.opto)
  fprintf('creating layout from cfg.opto\n');
  sens = cfg.opto;
  if (hasdata)
    layout = opto2lay(sens, data.label);
  else
    layout = opto2lay(sens, sens.label);
  end;
  
elseif isfield(data, 'opto') && isstruct(data.opto)
  fprintf('creating layout from data.opto\n');
  sens = data.opto;
  if (hasdata)
    layout = opto2lay(sens, data.label);
  else
    layout = opto2lay(sens, sens.label);
  end;
  
elseif (~isempty(cfg.image) || ~isempty(cfg.mesh)) && isempty(cfg.layout)
  % deal with image file
  if ~isempty(cfg.image)
    
    fprintf('reading background image from %s\n', cfg.image);
    [p,f,e] = fileparts(cfg.image);
    switch e
      case '.mat'
        tmp    = load(cfg.image);
        fnames = fieldnames(tmp);
        if numel(fnames)~=1
          error('there is not just a single variable in %s', cfg.image);
        else
          img = tmp.(fname{1});
        end
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
    
  elseif ~isempty(cfg.mesh)
    if isfield(cfg.mesh, 'sulc')
      ft_plot_mesh(cfg.mesh, 'edgecolor','none','vertexcolor',cfg.mesh.sulc);colormap gray;
    else
      ft_plot_mesh(cfg.mesh, 'edgecolor','none');
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
          if ~isempty(cfg.image)
            h = image(img);
          else
            h = ft_plot_mesh(cfg.mesh,'edgecolor','none','vertexcolor',cfg.mesh.sulc);
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
          if ~isempty(cfg.image)
            h = image(img);
          else
            h = ft_plot_mesh(cfg.mesh,'edgecolor','none','vertexcolor',cfg.mesh.sulc);
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
          if ~isempty(cfg.image)
            h = image(img);
          else
            h = ft_plot_mesh(cfg.mesh,'edgecolor','none','vertexcolor',cfg.mesh.sulc);
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

% set width and height if not already present

% FIXME there is a conflict between the use of cfg.style here and in topoplot
if ~strcmp(cfg.style, '3d')
  
  % input is 'scalp' data for which generating a head outline makes sense
  label = ft_channelselection(layout.label, setdiff(layout.label, {'COMNT', 'SCALE'}));
  if (ft_senstype(label,'meg') || ft_senstype(label,'eeg') || ft_senstype(label,'nris')) && (isempty(cfg.headshape) && isempty(cfg.mri))
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
  else % data is not 'scalp' level data and auto head outline+mask makes no sense. See if we can generate outline/mask based on something else
    
    % a mesh/mri was specified for outline generation (in case an outline already exist, add it to it)
    if (~isempty(cfg.headshape) || ~isempty(cfg.mri))
      % parse headshape/mri
      if ~isempty(cfg.headshape)
        if ischar(cfg.headshape) && exist(cfg.headshape, 'file')
          fprintf('reading headshape from file %s\n', cfg.headshape);
          outlbase = ft_read_headshape(cfg.headshape);
        elseif ft_datatype(cfg.headshape,'mesh')
          outlbase = cfg.headshape;
        else
          error('cfg.headshape does not contain headshape')
        end
        outlbase = ft_datatype_headmodel(outlbase);
      elseif ~isempty(cfg.mri)
        if ischar(cfg.mri) && exist(cfg.mri, 'file')
          fprintf('reading MRI from file %s\n', cfg.mri);
          outlbase = ft_read_mri(cfg.mri);
        elseif ft_datatype(cfg.mri,'volume')
          outlbase = cfg.mri;
        else
          error('cfg.mri does not contain mri')
        end
        outlbase = ft_datatype_volume(outlbase);
        % create mesh from anatomical field, and use as headshape below
        cfgpm = [];
        cfgpm.method      = 'projectmesh';
        cfgpm.tissue      = 'brain';
        cfgpm.numvertices = 1e5;
        outlbase = ft_prepare_mesh(cfgpm, outlbase);
      end
      
      % check for presence of cfg.viewpoint
      if isempty(cfg.viewpoint) 
        warning('cfg.headshape/mri is supplied without explicit orthographic projection viewpoint, assuming superior view')
        cfg.viewpoint = 'superior';
      end
      
      % check coordinate system of outlbase
      if isempty(outlbase.coordsys)
        error('no coordsys field found in mri/headshape, coordinate system unknown, use ft_determine_coordsys')
      end
      if exist('sens','var')
        if ~isfield(sens,'coordsys')
          error('no coordsys field found in elec/grad/opto, coordinate system unknown, use ft_determine_coordsys')
        elseif ~isequal(outlbase.coordsys,sens.coordsys)
          error('coordinate system of mri/headshape does not match that of sensors')
        end
      else
        % if used, sens would exist
        error('using cfg.headshape/mri to generate layout outline is not supported without providing elec/grad structure/file')
      end
      
      % match units with that of electrodes
      if exist('sens','var')
        outlbase = ft_convert_units(outlbase,sens.unit);
      end
      
      % generate outline based on matlab version
      if ft_platform_supports('boundary')
        
        % extract points indicating brain
        braincoords = outlbase.pos;
        % apply orthographic projection and extract XY
        braincoords = getorthoviewpos(braincoords,outlbase.coordsys,cfg.viewpoint);
        
        % get outline
        k = boundary(braincoords,.8);
        outline = braincoords(k,:);
        
      else % fallback, sad!
        
        % plot mesh in rotated view, rotate, screencap, and trace frame to generate outline
        outlbase.pos = getorthoviewpos(outlbase.pos,outlbase.coordsys,cfg.viewpoint);
        h = figure('visible','off');
        ft_plot_mesh(outlbase,'facecolor',[0 0 0], 'EdgeColor', 'none');
        view([0 90])
        axis tight
        xlim = get(gca,'xlim');
        ylim = get(gca,'ylim');
        set(gca,'OuterPosition',[0.2 0.2 .6 .6]) % circumvent weird matlab bug, wth?
        % extract frame for tracing
        drawnow % need to flush buffer, otherwise frame will not get extracted properly
        frame     = getframe(h);
        close(h)
        imtotrace = double(~logical(sum(frame.cdata,3))); % needs to be binary to trace
        % image is not in regular xy space, flip, and transpose
        imtotrace = flipud(imtotrace).';
        
        % trace image generated above
        [row col] = find(imtotrace,1,'first'); % set arbitrary starting point
        trace = bwtraceboundary(imtotrace,[row, col],'N');
        
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
        
        %
        outline = [x y];
      end
      
      % subsample outline
      if size(outline,1)>5e3 % 5e3 points should be more than enough to get an outine with acceptable detail
        outline = outline(1:floor(size(outline,1)/5e3):end,:);
      end
      
      % save outline
      if isfield(layout,'outline')
        layout.outline{end+1} = outline;
      else
        layout.outline = {outline};
      end
    end
    
    % Below, it is assumed that 'a mask' is always preferred to 'no mask'. In the worst case scenario, it would lead to ugly topoplots, if the
    % user requests them.
    % no mask is present, generate mask as boundary/convex hull of pos+width/height
    if ~isfield(layout, 'mask') || (isfield(layout,'mask') && isempty(layout.mask))
      % get index of non-COMNT/SCALE in case present
      ind1 = strcmp(layout.label,'COMNT');
      ind2 = strcmp(layout.label,'SCALE');
      ind  = ~(ind1 | ind2);
      % get position of all sensor 'boxes' and draw boundary around it, or, in older matlabs, a convex hull
      x = layout.pos(ind,1);
      y = layout.pos(ind,2);
      w = layout.width(ind);
      h = layout.height(ind);
      boxpos = [x - (w/2) y - (h/2);... % lb
        x - (w/2) y + (h/2);... % lt
        x + (w/2) y - (h/2);... % rb
        x + (w/2) y + (h/2)]; % rt
      if ft_platform_supports('boundary')
        k = boundary(boxpos,.5);
        mask = boxpos(k,:);
      else
        mask = boxpos(convhull(boxpos(:,1),boxpos(:,2)),:);
      end
      layout.mask{1} = mask;
    end
    
  end
end % if style=3d



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
  XY = cat(1, layout.pos, layout.outline{:}, layout.mask{:});
  layout.pos(end+1,:)  = [min(XY(:,1)) min(XY(:,2))];
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
  XY = cat(1, layout.pos, layout.outline{:}, layout.mask{:});
  layout.pos(end+1,:)  = [max(XY(:,1)) min(XY(:,2))];
elseif any(strcmp('SCALE', layout.label)) && skipscale
  % remove the scale entry
  sel = find(strcmp('SCALE', layout.label));
  layout.label(sel) = [];
  layout.pos(sel,:) = [];
  layout.width(sel) = [];
  layout.height(sel) = [];
end

% ensure proper format of some of label (see bug 1909 -roevdmei)
layout.label  = layout.label(:);

% to plot the layout for debugging, you can use this code snippet
if strcmp(cfg.feedback, 'yes') && strcmpi(cfg.style, '2d')
  tmpcfg = [];
  tmpcfg.layout = layout;
  ft_layoutplot(tmpcfg);
end

% to write the layout to a .mat or text file, you can use this code snippet
if ~isempty(cfg.output) && strcmpi(cfg.style, '2d')
  fprintf('writing layout to ''%s''\n', cfg.output);
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
  error('writing a 3D layout to an output file is not supported');
end


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
if isfield(sens, 'chanposorg')
  chanposorg = sens.chanposorg;
else
  chanposorg = [];
end
if isfield(sens, 'balance') && ~strcmp(sens.balance.current, 'none')
  sens = undobalancing(sens);
  if size(chanposorg, 1) == numel(sens.label)
    sens.chanpos = chanposorg;
  end
  % In case not all the locations have NaNs it might still be useful to plot them
  % But perhaps it'd be better to have any(any
elseif any(all(isnan(sens.chanpos)))
  [sel1, sel2] = match_str(sens.label, sens.labelorg);
  sens.chanpos = chanposorg(sel2, :);
  sens.label   = sens.labelorg(sel2);
end

fprintf('creating layout for %s system\n', ft_senstype(sens));

% apply rotation, but only if viewpoint is not used specifically
if isempty(viewpoint)
  if isempty(rotatez)
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
  sens.chanpos = ft_warp_apply(rotate([0 0 rotatez]), sens.chanpos, 'homogenous');
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
    % deal with coordinate system
    if isempty(sens.coordsys)
      error('no coordinate system found, elec/grad in unknown coordinate system, use ft_determine_coordsys')
    end
    
    % determine auto view if requested
    switch sens.coordsys
      case {'tal','mni'}
        if strcmp(viewpoint,'auto')
          % simple automatic determination of 'ideal' viewpoint
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
      otherwise
        error(['cfg.viewpoint = ''auto'' is not supported for coordinate system ' sens.coordsys])
    end
    
    % 3D to 2D
    prj = getorthoviewpos(pos,sens.coordsys,viewpoint);
  end
  
  % this copy will be used to determine the minimum distance between channels
  % we need a copy because prj retains the original positions, and
  % prjForDist might need to be changed if the user wants to keep
  % overlapping channels
  % also subselect channels for computing width/height if requested by cfg.boxchannel
  boxchannel = ft_channelselection(boxchannel, label); 
  boxchansel    = match_str(label, boxchannel); 
  prjForDist = prj(boxchansel,:);
  
  % check whether many channels occupy identical positions, if so shift
  % them around if requested
  if size(unique(prjForDist,'rows'),1) / size(prjForDist,1) < 0.8
    if strcmp(overlap, 'shift')
      ft_warning('the specified sensor configuration has many overlapping channels, creating a layout by shifting them around (use a template layout for better control over the positioning)');
      prj = shiftxy(prj', 0.2)';
      prjForDist = prj(boxchansel,:);
    elseif strcmp(overlap, 'no')
      error('the specified sensor configuration has many overlapping channels, you specified not to allow that');
    elseif strcmp(overlap, 'keep')
      prjForDist = unique(prj(boxchansel,:), 'rows');
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
  else
    mindist = eps; % not sure this is a good value but it's just to prevent crashes when
    % the EEG sensor definition is meaningless
  end
  
  X = prj(:,1);
  Y = prj(:,2);
  Width  = ones(size(X)) * mindist * 0.8;
  Height = ones(size(X)) * mindist * 0.6;
  layout.pos    = [X Y];
  layout.width  = Width;
  layout.height = Height;
  layout.pos    = prj;
  layout.label  = label;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% convert 2D optode positions into 2D layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = opto2lay(opto, label)
layout = [];
layout.pos = [];
layout.label = {};
layout.width = [];
layout.height = [];

[rxnames, rem] = strtok(label, {'-', ' '});
[txnames, rem] = strtok(rem, {'-', ' '});

for i=1:numel(label)
  % create average positions
  rxid = ismember(opto.fiberlabel, rxnames(i));
  txid = ismember(opto.fiberlabel, txnames(i));
  layout.pos(i, :) = opto.fiberpos(rxid, :)/2 + opto.fiberpos(txid, :)/2;
  layout.label(end+1)  = label(i);
  layout.width(end+1)  = 1;
  layout.height(end+1) = 1;
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% obtain XY pos from XYZ pos as orthographic projections dependent on viewpoint/coordsys 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getorthoviewpos(pos,coordsys,viewpoint)
% 
if size(pos,2)~=3
  error('XYZ coordinates needed to obtain orthographic projections based on viewpoint')
end

% create view(az,el) transformation matrix
switch coordsys
  case {'tal','mni'}
    switch viewpoint
      case 'left'
        transmat = viewmtx(-90,0); 
      case 'right'
        transmat = viewmtx(90,0); 
      case 'superior'
        transmat = viewmtx(0,90); 
      case 'inferior'
        transmat = viewmtx(180,-90); 
      case 'posterior'
        transmat = viewmtx(0,0); 
      case 'anterior'
        transmat = viewmtx(180,0); 
      otherwise
        error(['orthographic projection using viewpoint ' viewpoint ' is not supported'])
    end
  otherwise
    error(['orthographic projection using viewpoint ' viewpoint ' is not supported for coordinate system ' coordsys])
end

% extract xy
pos      = ft_warp_apply(transmat,pos,'homogenous');
pos      = pos(:,[1 2]);

