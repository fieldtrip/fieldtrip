function [elec] = ft_electrodeplacement(cfg, varargin)

% FT_ELECTRODEPLACEMENT allows manual placement of electrodes on a MRI scan, CT scan
% or on a triangulated surface of the head. This function supports different methods.
%
% VOLUME - Navigate an orthographic display of a volume (e.g. CT or MRI scan), and
% assign an electrode label to the current crosshair location by clicking on a label
% in the eletrode list. You can undo the selection by clicking on the same label
% again. The electrode labels shown in the list can be prespecified using cfg.channel
% when calling ft_electrodeplacement. The zoom slider allows zooming in at the
% location of the crosshair. The intensity sliders allow thresholding the image's low
% and high values. The magnet feature transports the crosshair to the nearest peak
% intensity voxel, within a certain voxel radius of the selected location. The labels
% feature displays the labels of the selected electrodes within the orthoplot. The
% global feature allows toggling the view between all and near-crosshair
% markers. The scan feature allows toggling between scans when another scan
% is given as input.
%
% HEADSHAPE - Navigate a triangulated scalp (for EEG) or brain (for ECoG) surface,
% and assign an electrode location by clicking on the surface. The electrode is
% placed on the triangulation itself.
%
% 1020 - Starting from a triangulated scalp surface and the nasion, inion, left and
% right pre-auricular points, this automatically constructs and follows contours over
% the surface according to the 5% system. Electrodes are placed at certain relative
% distances along these countours. This is an extension of the 10-20 standard
% electrode placement system and includes the 20%, 10% and 5% locations. See
% "Oostenveld R, Praamstra P. The five percent electrode system for high-resolution
% EEG and ERP measurements. Clin Neurophysiol. 2001 Apr;112(4):713-9" for details.
%
% SHAFT - This is for placing electrodes along a linear sEEG shaft. The tip of the
% shaft corresponding to the first electrode, another point along the shaft, and the
% distance between the electrodes should be specified. If the shaft is not straight
% but curved, you should specify multiple (at least two) points along the shaft,
% i.e., specify cfg.shaft.along as an Nx3 array for N points along the shaft. The
% number of electrodes to be distributed along the shaft is determined from cfg.channel.
%
% GRID - This is for placing electrodes on a regular MxN ECoG grid. Each of the four
% cornerpoints of the grid must be specified, along with the dimensions of the grid.
% Following piecewise linear placement of the electrodes on the grid, you can use
% FT_ELECTRODEREALIGN with cfg.method='project' to project them to the curved brain
% surface.
%
% Use as
%   [elec] = ft_electrodeplacement(cfg, mri)
%   [elec] = ft_electrodeplacement(cfg, ct)
%   [elec] = ft_electrodeplacement(cfg, mri, ct, ..)
% where the second and subsequent input arguments should be one or multiple
% anatomical MRIs and/or CTs, or
%   [elec] = ft_electrodeplacement(cfg, headshape)
% where the input headshape should be a surface triangulation.
%
% The configuration can contain the following options
%   cfg.method         = string representing the method for placing the electrodes
%                        'volume'          interactively locate electrodes on three orthogonal slices of a volumetric MRI or CT scan
%                        'headshape'       interactively locate electrodes on a head surface
%                        '1020'            automatically locate electrodes on a head surface according to the 10-20 system
%                        'shaft'           automatically locate electrodes along a linear sEEG shaft
%                        'grid'            automatically locate electrodes on a MxN ECoG grid
%   cfg.figure         = 'yes' or 'no', whether to open a new figure. You can also specify a figure handle from FIGURE, GCF or SUBPLOT. (default = 'yes')
%   cfg.position       = location and size of the figure, specified as [left bottom width height] (default is automatic)
%   cfg.renderer       = string, 'opengl', 'zbuffer', 'painters', see RENDERERINFO (default = 'opengl')
%
% The following options apply to the 'volume' method
%   cfg.parameter      = string, field in data (default = 'anatomy' if present in data)
%   cfg.channel        = Nx1 cell-array with selection of channels (default = {'1' '2' ...})
%   cfg.elec           = struct containing previously placed electrodes (this overwrites cfg.channel)
%   cfg.clim           = color range of the data (default = [0 1], i.e. the full range)
%   cfg.magtype        = string representing the 'magnet' type used for placing the electrodes
%                        'peakweighted'    place electrodes at weighted peak intensity voxel (default)
%                        'troughweighted'  place electrodes at weighted trough intensity voxel
%                        'peak'            place electrodes at peak intensity voxel (default)
%                        'trough'          place electrodes at trough intensity voxel
%                        'weighted'        place electrodes at center-of-mass
%   cfg.magradius      = number representing the radius for the cfg.magtype based search (default = 3)
%
% The following options apply to the '1020' method
%   cfg.fiducial.nas   = 1x3 vector with coordinates
%   cfg.fiducial.ini   = 1x3 vector with coordinates
%   cfg.fiducial.lpa   = 1x3 vector with coordinates
%   cfg.fiducial.rpa   = 1x3 vector with coordinates
%   cfg.feedback       = string, can be 'yes' or 'no' for detailed feedback (default = 'yes')
%
% The following options apply to the 'shaft' method
%   cfg.shaft.tip      = 1x3 position of the electrode at the tip of the shaft
%   cfg.shaft.along    = 1x3 or Nx3 positions along the shaft
%   cfg.shaft.distance = scalar, distance between electrodes
%
% The following options apply to the 'grid' method
%   cfg.grid.corner1   = 1x3 position of the upper left corner point
%   cfg.grid.corner2   = 1x3 position of the upper right corner point
%   cfg.grid.corner3   = 1x3 position of the lower left corner point
%   cfg.grid.corner4   = 1x3 position of the lower right corner point
%
% In the interactive 'headshape' and 'volume' methods you can click once on an
% already assigned electrode to jump to that electrode position and you can click
% twice to remove the assigned electrode position.
%
% See also FT_ELECTRODEREALIGN, FT_VOLUMEREALIGN, FT_VOLUMESEGMENT, FT_PREPARE_MESH

% Copyright (C) 2015-2019, Arjen Stolk, Sandon Griffin & Robert Oostenveld
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
ft_preamble loadvar mri
ft_preamble provenance mri

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% ensure that old and unsupported options are not being relied on by the end-user's script
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2837

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamed',    {'viewdim', 'axisratio'});
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'mri', 'volume'});
cfg = ft_checkconfig(cfg, 'renamed',    {'newfigure', 'figure'});

% set the defaults
cfg.method        = ft_getopt(cfg, 'method',              []); % volume, headshape, 1020, shaft
cfg.feedback      = ft_getopt(cfg, 'feedback',         'yes');
cfg.parameter     = ft_getopt(cfg, 'parameter',    'anatomy');
cfg.channel       = ft_getopt(cfg, 'channel',             []); % default will be determined further down {'1', '2', ...}
cfg.elec          = ft_getopt(cfg, 'elec',                []); % use previously placed electrodes
cfg.flip          = ft_getopt(cfg, 'flip',                []); % the default is set below
cfg.renderer      = ft_getopt(cfg, 'renderer',      'opengl');
% view options
cfg.clim          = ft_getopt(cfg, 'clim',             [0 1]); % initial volume intensity limit voxels
cfg.markerdist    = ft_getopt(cfg, 'markerdist',           5); % marker-slice distance view when ~global
% magnet options
cfg.magtype       = ft_getopt(cfg, 'magtype', 'peakweighted'); % detect weighted peaks or troughs
cfg.magradius     = ft_getopt(cfg, 'magradius',            3); % specify the physical unit radius
cfg.voxelratio    = ft_getopt(cfg, 'voxelratio',      'data'); % display size of the voxel, 'data' or 'square'
cfg.axisratio     = ft_getopt(cfg, 'axisratio',       'data'); % size of the axes of the three orthoplots, 'square', 'voxel', or 'data'

if isempty(cfg.method) && ~isempty(varargin)
  % the default determines on the input data
  switch ft_datatype(varargin{1})
    case 'volume'
      cfg.method = 'volume';
    case 'mesh'
      cfg.method = 'headshape';
    case 'source+mesh'
      cfg.method = 'headshape';
  end
end

if isempty(cfg.flip)
  if strcmp(cfg.method, 'volume')
    cfg.flip = 'yes';
  else
    cfg.flip = 'no';
  end
end

% check if the input data is valid for this function
switch cfg.method
  case 'volume'
    for v = 1:numel(varargin)
      mri{v} = ft_checkdata(varargin{v}, 'datatype', 'volume', 'feedback', 'yes', 'hascoordsys', 'yes', 'hasunit', 'yes');
    end
  case  {'headshape', '1020'}
    headshape = fixpos(varargin{1});
    headshape = ft_checkdata(headshape, 'hascoordsys', 'yes', 'hasunit', 'yes');
end

if strcmp(cfg.flip, 'yes')
  for v = 1:numel(varargin)
    % align the anatomical volume approximately to coordinate system, this puts it upright
    origmethod = cfg.method;
    tmpcfg = [];
    tmpcfg.method = 'flip';
    tmpcfg.trackcallinfo = 'no';
    tmpcfg.showcallinfo = 'no';
    mri{v} = ft_volumereslice(tmpcfg, mri{v});
    [cfg, mri] = rollback_provenance(cfg, mri);
    cfg.method = origmethod;
  end
end

% setup electrode labels if possible
chanlabel = {}; 
chanstring = {};
markerlab = {}; 
markerpos = {};
if ~isempty(cfg.elec) % re-use previously placed (cfg.elec) electrodes
  for e = 1:numel(cfg.elec.label)
    chanlabel{end+1,1} = cfg.elec.label{e};
    chanstring{end+1} = ['<HTML><FONT color="black">' cfg.elec.label{e} '</FONT></HTML>']; % hmtl'ize

    markerlab{end+1,1} = cfg.elec.label{e};
    markerpos{end+1,1} = cfg.elec.elecpos(e,:);
  end
end
if ~isempty(cfg.channel) % use prespecified (cfg.channel) electrode labels
  for c = 1:numel(cfg.channel)
    if ~ismember(cfg.channel{c}, chanlabel) % avoid overlap between cfg.channel and elec.label
      chanlabel{end+1,1} = cfg.channel{c};
      chanstring{end+1} = ['<HTML><FONT color="silver">' cfg.channel{c} '</FONT></HTML>']; % hmtl'ize

      markerlab{end+1,1} = [];
      markerpos{end+1,1} = zeros(0,3);
    end
  end
end
if isempty(cfg.elec) && isempty(cfg.channel) % create electrode labels on-the-fly
  for c = 1:300
    chanlabel{end+1,1} = sprintf('%d', c);
    chanstring{end+1} = ['<HTML><FONT color="silver">' sprintf('%d', c) '</FONT></HTML>']; % hmtl'ize

    markerlab{end+1,1} = {};
    markerpos{end+1,1} = zeros(0,3);
  end
end

% this is where the different methods are implemented
switch cfg.method

  case 'volume'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is an interactive method, start building the figure

    h = open_figure(keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'}));
    set(h, 'Name', mfilename);
    set(h, 'Units', 'normalized');
    set(h, 'Color', [1 1 1]);
    set(h, 'MenuBar', 'none');
    set(h, 'windowbuttondownfcn', @cb_buttonpress);
    set(h, 'windowbuttonupfcn',   @cb_buttonrelease);
    set(h, 'windowkeypressfcn',   @cb_keyboard);
    set(h, 'CloseRequestFcn',     @cb_quit);

    % volume-dependent axis settings
    for v = 1:numel(mri)
      if strcmp(cfg.axisratio, 'voxel')
        % determine the number of voxels to be plotted along each axis
        axlen1 = mri{v}.dim(1);
        axlen2 = mri{v}.dim(2);
        axlen3 = mri{v}.dim(3);
      elseif strcmp(cfg.axisratio, 'data')
        % determine the length of the edges along each axis
        [cp_voxel, cp_head] = cornerpoints(mri{v}.dim, mri{v}.transform);
        axlen1 = norm(cp_head(2,:)-cp_head(1,:));
        axlen2 = norm(cp_head(4,:)-cp_head(1,:));
        axlen3 = norm(cp_head(5,:)-cp_head(1,:));
      elseif strcmp(cfg.axisratio, 'square')
        % the length of the axes should be equal
        axlen1 = 1;
        axlen2 = 1;
        axlen3 = 1;
      end

      % this is the size reserved for subplot h1, h2 and h3
      h1size(1) = 0.92*axlen1/(axlen1 + axlen2); % x
      h1size(2) = 0.92*axlen3/(axlen2 + axlen3); % z
      h2size(1) = 0.92*axlen2/(axlen1 + axlen2); % y
      h2size(2) = 0.92*axlen3/(axlen2 + axlen3); % z
      h3size(1) = 0.92*axlen1/(axlen1 + axlen2); % x
      h3size(2) = 0.92*axlen2/(axlen2 + axlen3); % y

      % axis handles will hold the anatomical functional if present, along with labels etc.
      h1 = axes('position', [0.02                0.02+0.04+h3size(2) h1size(1) h1size(2)]); % x z
      h2 = axes('position', [0.02+0.04+h1size(1) 0.02+0.04+h3size(2) h2size(1) h2size(2)]); % y z
      h3 = axes('position', [0.02                0.02                h3size(1) h3size(2)]); % x y

      set(h1, 'Tag', 'ik', 'Visible', 'off', 'XAxisLocation', 'top'); axis(h1, 'equal');
      set(h2, 'Tag', 'jk', 'Visible', 'off', 'YAxisLocation', 'right'); axis(h2, 'equal'); % after rotating in ft_plot_ortho this becomes top
      set(h3, 'Tag', 'ij', 'Visible', 'off'); axis(h3, 'equal');

      if strcmp(cfg.voxelratio, 'square')
        voxlen1 = 1;
        voxlen2 = 1;
        voxlen3 = 1;
      elseif strcmp(cfg.voxelratio, 'data')
        % the size of the voxel is scaled with the data
        [cp_voxel, cp_head] = cornerpoints(mri{v}.dim, mri{v}.transform);
        voxlen1 = norm(cp_head(2,:)-cp_head(1,:))/norm(cp_voxel(2,:)-cp_voxel(1,:));
        voxlen2 = norm(cp_head(4,:)-cp_head(1,:))/norm(cp_voxel(4,:)-cp_voxel(1,:));
        voxlen3 = norm(cp_head(5,:)-cp_head(1,:))/norm(cp_voxel(5,:)-cp_voxel(1,:));
      end

      %set(h1, 'DataAspectRatio', 1./[voxlen1 voxlen2 voxlen3]); % FIXME: this no longer works when using mri.transform with ft_plot_ortho (instead of eye(4));
      %set(h2, 'DataAspectRatio', 1./[voxlen1 voxlen2 voxlen3]);
      %set(h3, 'DataAspectRatio', 1./[voxlen1 voxlen2 voxlen3]);

      mri{v}.axes = [h1 h2 h3];
      mri{v}.h1size = h1size;
      mri{v}.h2size = h2size;
      mri{v}.h3size = h3size;
      mri{v}.clim = cfg.clim;
      mri{v}.slim = [.9 1]; % 90% - maximum

      dat = double(mri{v}.(cfg.parameter));
      dmin = min(dat(:));
      dmax = max(dat(:));
      mri{v}.dat = (dat-dmin)./(dmax-dmin); % range between 0 and 1
      clear dat dmin dmax
    end

    % intensity range sliders
    h45text = uicontrol('Style', 'text',...
      'String', 'Intensity',...
      'Units', 'normalized', ...
      'Position', [2*mri{1}.h1size(1)-0.02 mri{1}.h3size(2)-0.02 mri{1}.h1size(1)/4 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on');

    h4 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 1, ...
      'Value', cfg.clim(1), ...
      'Units', 'normalized', ...
      'Position', [2*mri{1}.h1size(1)-0.03 0.10+mri{1}.h3size(2)/3 0.05 mri{1}.h3size(2)/2-0.05], ...
      'Callback', @cb_minslider);

    h5 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 1, ...
      'Value', cfg.clim(2), ...
      'Units', 'normalized', ...
      'Position', [2*mri{1}.h1size(1)+0.02 0.10+mri{1}.h3size(2)/3 0.05 mri{1}.h3size(2)/2-0.05], ...
      'Callback', @cb_maxslider);

    % java intensity range slider (dual-knob slider): the java component gives issues when wanting to
    % access the opt structure
    % [jRangeSlider] = com.jidesoft.swing.RangeSlider(0,1,cfg.clim(1),cfg.clim(2));  % min,max,low,high
    % [jRangeSlider, h4] = javacomponent(jRangeSlider, [], h);
    % set(h4, 'Units', 'normalized', 'Position', [0.05+h1size(1) 0.07 0.07 h3size(2)], 'Parent', h);
    % set(jRangeSlider, 'Orientation', 1, 'PaintTicks', true, 'PaintLabels', true, ...
    %     'Background', java.awt.Color.white, 'StateChangedCallback', @cb_intensityslider);

    % electrode listbox
    h6 = uicontrol('Style', 'listbox', ...
      'Parent', h, ...
      'Value', [], 'Min', 0, 'Max', numel(chanstring), ...
      'Units', 'normalized', ...
      'FontSize', 12, ...
      'Position', [mri{1}.h1size(1)+0.07 0.02 mri{1}.h2size(1)/2.5 mri{1}.h3size(2)], ...
      'Callback', @cb_eleclistbox, ...
      'String', chanstring);

    % switches / radio buttons
    h7text = uicontrol('Style', 'text',...
      'String', 'Magnet',...
      'Units', 'normalized', ...
      'Position', [2*mri{1}.h1size(1)-0.047 0.18 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on');

    h7 = uicontrol('Style', 'popupmenu',...
      'Parent', h, ...
      'Value', 4, ... % corresponding to magradius = 3 (see String)
      'String', {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}, ...
      'Units', 'normalized', ...
      'Position', [2*mri{1}.h1size(1)-0.103 0.18 mri{1}.h1size(1)/4.25 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Callback', @cb_magnetbutton);
    radii = get(h7, 'String');
    if ~ismember(num2str(cfg.magradius), radii) % add user-specified radius to the list
      set(h7, 'String', [radii(:); num2str(cfg.magradius)]);
      set(h7, 'Value', numel(radii)+1);
    end

    h8 = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 0, ...
      'String', 'Labels',...
      'Units', 'normalized', ...
      'Position', [2*mri{1}.h1size(1)-0.05 0.14 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Callback', @cb_labelsbutton);

    h9 = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 0, ...
      'String', 'Global',...
      'Units', 'normalized', ...
      'Position', [2*mri{1}.h1size(1)-0.05 0.10 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Callback', @cb_globalbutton);

    h11 = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 0, ...
      'String', 'Scatter',...
      'Units', 'normalized', ...
      'Position', [2*mri{1}.h1size(1)-0.05 0.06 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Callback', @cb_scatterbutton);

    h12 = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 0, ...
      'String', 'CT/MRI',...
      'Units', 'normalized', ...
      'Position', [2*mri{1}.h1size(1)-0.05 0.02 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Visible', 'off', ...
      'Callback', @cb_scanbutton);
    if numel(mri)>1; set(h12, 'Visible', 'on'); end % only when two scans are given as input

    % zoom slider
    h10text = uicontrol('Style', 'text',...
      'String', 'Zoom',...
      'Units', 'normalized', ...
      'Position', [1.8*mri{1}.h1size(1)-0.04 mri{1}.h3size(2)-0.02 mri{1}.h1size(1)/4 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on');

    h10 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 0.9, ...
      'Value', 0, ...
      'Units', 'normalized', ...
      'Position', [1.8*mri{1}.h1size(1)-0.03 0.10+mri{1}.h3size(2)/3 0.05 mri{1}.h3size(2)/2-0.05], ...
      'SliderStep', [.1 .1], ...
      'Callback', @cb_zoomslider);

    % create structure to be passed to gui
    opt               = [];
    opt.method        = 'volume'; % this is to distinguish between volume and headshape in the callbacks
    opt.label         = chanlabel;
    opt.handlesaxes   = [mri{1}.axes(1) mri{1}.axes(2) mri{1}.axes(3) h4 h5 h6 h7 h8 h9 h10 h11 h12];
    opt.mainfig       = h;
    opt.quit          = false;
    opt.update        = [1 1 1];
    opt.init          = true;
    opt.tag           = 'ik';
    opt.ana           = mri{1}.dat; % the plotted anatomy
    opt.mri           = mri;
    opt.currmri       = 1;
    opt.showcrosshair = true;
    opt.pos           = [0 0 0]; % middle of the scan, head coordinates
    opt.showsurface   = true;
    opt.showlabels    = false;
    opt.showcolors    = false;
    opt.magnet        = get(h7, 'Value');
    opt.magradius     = cfg.magradius;
    opt.magtype       = cfg.magtype;
    opt.showmarkers   = true;
    opt.global        = get(h9, 'Value'); % show all markers in the current slices
    opt.scatter       = get(h11, 'Value'); % additional scatterplot
    opt.scan          = get(h12, 'Value'); % switch scans
    opt.slim          = [.9 1]; % 90% - maximum
    opt.markerlab     = markerlab;
    opt.markerpos     = markerpos;
    opt.markerdist    = cfg.markerdist; % hidden option
    opt.clim          = cfg.clim;
    opt.zoom          = 0;

    setappdata(h, 'opt', opt);
    cb_redraw(h);
    cb_help(h);

    while(opt.quit==0)
      uiwait(h);
      opt = getappdata(h, 'opt');
    end
    delete(h);

    % collect the results
    elec = keepfields(mri{1}, {'unit', 'coordsys'});
    elec.label    = {};
    elec.elecpos  = [];
    for i=1:length(opt.markerlab)
      if ~isempty(opt.markerlab{i,1})
        elec.label = [elec.label; opt.markerlab{i,1}];
        elec.elecpos = [elec.elecpos; opt.markerpos{i,1}];
      end
    end
    elec.chanpos = elec.elecpos;
    elec.tra = eye(size(elec.elecpos,1));

  case 'headshape'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is an interactive method, start building the figure

    h = open_figure(keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'}));
    set(h, 'Name', mfilename);
    set(h, 'Units', 'normalized');
    set(h, 'Color', [1 1 1]);
    set(h, 'windowbuttondownfcn', @cb_buttonpress);
    set(h, 'windowbuttonupfcn',   @cb_buttonrelease);
    set(h, 'windowkeypressfcn',   @cb_keyboard);
    set(h, 'CloseRequestFcn',     @cb_quit);

    % add an axis and rotate to the starting viewpoint
    h1 = axes(h, 'Visible', 'off'); view([90, 0]); axis off

    % electrode listbox
    h2 = uicontrol('Style', 'listbox', ...
      'Parent', h, ...
      'Value', [], 'Min', 0, 'Max', numel(chanstring), ...
      'Units', 'normalized', ...
      'FontSize', 12, ...
      'Position', [.8 0.001 .2 .5], ...
      'Callback', @cb_eleclistbox, ...
      'String', chanstring);

    h7 = uicontrol('Style', 'checkbox', ...
      'Parent', h, ...
      'Value', 1, ...
      'String', 'Surface', ...
      'Units', 'normalized', ...
      'Position', [.8 0.75 .2 .05], ...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Callback', @cb_surfacebutton);

    h8 = uicontrol('Style', 'checkbox', ...
      'Parent', h, ...
      'Value', 0, ...
      'String', 'Labels', ...
      'Units', 'normalized', ...
      'Position', [.8 0.70 .2 .05], ...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Callback', @cb_labelsbutton);

    h9 = uicontrol('Style', 'checkbox', ...
      'Parent', h, ...
      'Value', 1, ...
      'String', 'Colors', ...
      'Units', 'normalized', ...
      'Position', [.8 0.65 .2 .05], ...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Callback', @cb_colorsbutton);

    if ~isfield(headshape, 'color')
      set(h9, 'Visible', false);
    end

    h10 = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 1, ...
      'String', 'Camlight', ...
      'Units', 'normalized', ...
      'Position', [.8 0.60 .2 .05],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Callback', @cb_camlight);

    h11 = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 1, ...
      'String', 'Shiny', ...
      'Units', 'normalized', ...
      'Position', [.8 0.55 .2 .05],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility', 'on', ...
      'Callback', @cb_shiny);

    % create structure to be passed to gui
    opt               = [];
    opt.method        = 'headshape'; % this is to distinguish between volume and headshape in the callbacks
    opt.headshape     = headshape;
    opt.label         = chanlabel;
    opt.mainfig       = h;
    opt.quit          = false;
    opt.init          = true;
    opt.pos           = [0 0 0]; % middle of the scan, head coordinates (FIXME: this might mess up vertex finding, being an anchor)
    opt.showcrosshair = true;
    opt.showsurface   = true;
    opt.showlabels    = false;
    opt.showcolors    = true;
    opt.showmarkers   = true;
    opt.markerlab     = markerlab;
    opt.markerpos     = markerpos;
    opt.markerdist    = cfg.markerdist; % hidden option
    opt.camlight      = get(h10, 'value');
    opt.shiny         = get(h11, 'value');

    setappdata(h, 'opt', opt);
    cb_help(h);
    cb_redraw(h);
    view([90, 0]);

    while(opt.quit==0)
      uiwait(h);
      opt = getappdata(h, 'opt');
    end
    delete(h);

    % collect the results
    elec = keepfields(headshape, {'unit', 'coordsys'});
    elec.label    = {};
    elec.elecpos  = [];
    for i=1:length(opt.markerlab)
      if ~isempty(opt.markerlab{i,1})
        elec.label = [elec.label; opt.markerlab{i,1}];
        elec.elecpos = [elec.elecpos; opt.markerpos{i,1}];
      end
    end
    elec.chanpos = elec.elecpos;
    elec.tra = eye(size(elec.elecpos,1));

  case '1020'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is an automatic method without figure

    % the placement procedure fails if the fiducials coincide with vertices
    dist = @(x, y) sqrt(sum(bsxfun(@minus, x, y).^2,2));
    tolerance = 0.1 * ft_scalingfactor('mm', headshape.unit);  % 0.1 mm
    nas = cfg.fiducial.nas;
    ini = cfg.fiducial.ini;
    lpa = cfg.fiducial.lpa;
    rpa = cfg.fiducial.rpa;
    if any(dist(headshape.pos, nas)<tolerance)
      ft_warning('Nasion coincides with headshape vertex, adding random displacement of about %f %s', tolerance, headshape.unit);
      nas = nas + tolerance*randn(1,3);
    end
    if any(dist(headshape.pos, ini)<tolerance)
      ft_warning('Inion coincides with headshape vertex, adding random displacement of about %f %s', tolerance, headshape.unit);
      ini = ini + tolerance*randn(1,3);
    end
    if any(dist(headshape.pos, lpa)<tolerance)
      ft_warning('LPA coincides with headshape vertex, adding random displacement of about %f %s', tolerance, headshape.unit);
      lpa = lpa + tolerance*randn(1,3);
    end
    if any(dist(headshape.pos, rpa)<tolerance)
      ft_warning('RPA coincides with headshape vertex, adding random displacement of about %f %s', tolerance, headshape.unit);
      rpa = rpa + tolerance*randn(1,3);
    end

    % place the electrodes automatically according to the fiducials
    [pos, lab] = elec1020_locate(headshape.pos, headshape.tri, nas, ini, lpa, rpa, istrue(cfg.feedback));
    % construct the output
    elec = keepfields(headshape, {'unit', 'coordsys'});
    elec.elecpos = pos;
    elec.label   = lab(:);

  case 'shaft'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is an automatic method without figure

    if size(cfg.shaft.along,1)==1
      % two points are specified, use a simple linear method
      pos = cfg.shaft.tip;
      ori = cfg.shaft.along - cfg.shaft.tip;
      ori = ori/norm(ori);

      elec = [];
      elec.label = cfg.channel;
      for i=1:numel(elec.label)
        elec.elecpos(i,:) = pos + (i-1) * ori * cfg.shaft.distance;
      end

    else
      % more than two points are specified, use a spline
      pos = [
        cfg.shaft.tip
        cfg.shaft.along
        ];
      dist = sqrt(sum(diff(pos,1,1).^2, 2)); % distance between subsequent positions
      dist = [0; cumsum(dist)];              % cumulative distance from the tip
      cv = csapi(dist', pos');               % spline that goes through all positions

      elec = [];
      elec.label = cfg.channel;
      for i=1:numel(elec.label)
        dist = (i-1) * cfg.shaft.distance;
        elec.elecpos(i,:) = cfg.shaft.tip + fnval(cv, dist)';
      end
    end

  case 'grid'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is an automatic method without figure

    % The four cornerpoints specified by the user will not exactly be a rectangle,
    % but will have in plane distortion and out-of-plane twist.

    % Each cornerpoint with the two edges to its direct neighbours describes a parallellogram.
    % We compute the position of each electrode according to this parallellogram
    % and then compute the weighted average based on distance to the cornerpoints.

    % The cornerpoints are specified like this
    %
    %     1     2
    %     3     4

    % parallellogram starting in corner 1
    dir12 = cfg.grid.corner2-cfg.grid.corner1; dir12 = dir12/norm(dir12);
    dir13 = cfg.grid.corner3-cfg.grid.corner1; dir13 = dir13/norm(dir13);

    % parallellogram starting in corner 2
    dir21 = cfg.grid.corner1-cfg.grid.corner2; dir21 = dir21/norm(dir21);
    dir24 = cfg.grid.corner4-cfg.grid.corner2; dir24 = dir24/norm(dir24);

    % parallellogram starting in corner 3
    dir31 = cfg.grid.corner1-cfg.grid.corner3; dir31 = dir31/norm(dir31);
    dir34 = cfg.grid.corner4-cfg.grid.corner3; dir34 = dir34/norm(dir34);

    % parallellogram starting in corner 4
    dir42 = cfg.grid.corner2-cfg.grid.corner4; dir42 = dir42/norm(dir42);
    dir43 = cfg.grid.corner3-cfg.grid.corner4; dir43 = dir43/norm(dir43);

    % assuming a 3x4 grid, the electrodes will be placed like this
    %
    %    1  2  3  4
    %    5  6  7  8
    %    9 10 11 12

    % specify the number of electrodes in the vertical and horizontal direction
    nv = cfg.grid.dim(1);
    nh = cfg.grid.dim(2);

    % determine the inter-electrode distance in the vertical and horizontal direction
    dv = (norm(cfg.grid.corner3-cfg.grid.corner1)/(nv-1) + norm(cfg.grid.corner4-cfg.grid.corner2)/(nv-1))/2;
    dh = (norm(cfg.grid.corner2-cfg.grid.corner1)/(nh-1) + norm(cfg.grid.corner4-cfg.grid.corner3)/(nh-1))/2;

    ft_notice('electrode spacing along 1st dimension is %f', dv);
    ft_notice('electrode spacing along 2nd dimension is %f', dh);

    % scale the
    dir12 = dir12*dh;
    dir13 = dir13*dv;
    dir21 = dir21*dh;
    dir24 = dir24*dv;
    dir31 = dir31*dv;
    dir34 = dir34*dh;
    dir42 = dir42*dv;
    dir43 = dir43*dh;

    if isempty(cfg.channel)
      ft_notice('using automatic channel labels');
      cfg.channel = arrayfun(@num2str, 1:(nh*nv), 'UniformOutput', false);
    else
      assert(numel(cfg.channel)==nh*nv, 'mismatch between cfg.grid.dim and cfg.channel');
    end

    pos1 = nan(nh*nv,3);
    pos2 = nan(nh*nv,3);
    pos3 = nan(nh*nv,3);
    pos4 = nan(nh*nv,3);

    % determine the position of each electrode according to each parallellogram
    k = 1;
    for i=1:nv
      for j=1:nh
        pos1(k,:) = ( i-1)*dir13 + ( j-1)*dir12 + cfg.grid.corner1;
        pos2(k,:) = ( i-1)*dir24 + (nh-j)*dir21 + cfg.grid.corner2;
        pos3(k,:) = (nv-i)*dir31 + ( j-1)*dir34 + cfg.grid.corner3;
        pos4(k,:) = (nv-i)*dir42 + (nh-j)*dir43 + cfg.grid.corner4;
        k = k+1;
      end % horizontal
    end % vertical

    d1 = nan(nh*nv,1);
    d2 = nan(nh*nv,1);
    d3 = nan(nh*nv,1);
    d4 = nan(nh*nv,1);

    % determine the schematic distance (i.e. counted in steps) to each corner
    k = 1;
    for i=1:nv
      for j=1:nh
        d1(k) = norm([( i-1)*dv ( j-1)*dh]);
        d2(k) = norm([( i-1)*dv (nh-j)*dh]);
        d3(k) = norm([(nv-i)*dv ( j-1)*dh]);
        d4(k) = norm([(nv-i)*dv (nh-j)*dh]);
        k = k+1;
      end % horizontal
    end % vertical

    % the weight is relative to the range, which is 0.5 times the dimagonal distance
    range = sqrt(((nv-1)*dv)^2 + ((nh-1)*dh)^2) / 2;
    % compute the weights from the distance
    d1 = 10.^(-d1/range);
    d2 = 10.^(-d2/range);
    d3 = 10.^(-d3/range);
    d4 = 10.^(-d4/range);
    % normalize the weights
    dd = (d1 + d2 + d3 + d4);
    d1 = d1 ./ dd;
    d2 = d2 ./ dd;
    d3 = d3 ./ dd;
    d4 = d4 ./ dd;

    % compute the weighted position
    elec = [];
    elec.label = cfg.channel;
    elec.elecpos = [d1 d1 d1].*pos1 + [d2 d2 d2].*pos2 + [d3 d3 d3].*pos3 + [d4 d4 d4].*pos4;

  otherwise
    ft_error('unsupported method ''%s''', cfg.method);

end % switch method

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   mri
ft_postamble provenance elec
ft_postamble history    elec
ft_postamble savevar    elec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_help(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

switch opt.method
  case 'volume'
    disp('==================================================================================');
    disp('0. Press "h" to show this help');
    disp('1. Options for viewing:');
    disp('   a. use the left mouse button to navigate in the image, or');
    disp('   b. use the arrow keys to increase or decrease the slice number by one');
    disp('2. Options for electrode placement:');
    disp('   a. click an electrode label in the list to assign it to the crosshair location, or');
    disp('   b. doubleclick a previously assigned electrode label to remove its marker');
    disp('3. Press "q" on the keyboard or close the window when you are done');
    disp('4. See Stolk, Griffin et al. Nature Protocols 2018 for further electrode processing options');

  case 'headshape'
    disp('==================================================================================');
    disp('0. Press "h" to show this help');
    disp('1. Options for electrode placement:');
    disp('   a. Use the mouse to click on the desired position for an electrode');
    disp('   b. Click on the corresponding electrode label');
    disp('   c. Double-click on an electrode label to remove the marker');
    disp('2. Press "v" to update the light position');
    disp('3. Press "q" or close the window when you are done');

end % switch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

if nargin<2
  eventdata = [];
end

switch opt.method
  case 'volume'
    cb_redraw_volume(h, eventdata)
  case 'headshape'
    cb_redraw_headshape(h, eventdata)
end % switch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw_volume(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

% determine current axes
set(0, 'CurrentFigure', opt.mainfig);
curr_ax = get(h, 'currentaxes');
tag = get(curr_ax, 'tag');

h1 = opt.handlesaxes(1);
h2 = opt.handlesaxes(2);
h3 = opt.handlesaxes(3);

if opt.init
  delete(findobj(opt.mainfig, 'Type', 'Surface')); % get rid of old orthos (to facilitate switching scans)
  ft_plot_ortho(opt.ana, 'transform', opt.mri{opt.currmri}.transform, 'coordsys', opt.mri{opt.currmri}.coordsys, 'location', opt.pos, 'style', 'subplot', 'parents', [h1 h2 h3], 'update', opt.update, 'doscale', false, 'clim', opt.clim, 'unit', opt.mri{opt.currmri}.unit);

  opt.anahandles = findobj(opt.mainfig, 'Type', 'Surface')';
  parenttag = get(opt.anahandles, 'parent');
  parenttag{1} = get(parenttag{1}, 'tag');
  parenttag{2} = get(parenttag{2}, 'tag');
  parenttag{3} = get(parenttag{3}, 'tag');
  [i1,i2,i3] = intersect(parenttag, {'ik';'jk';'ij'});
  opt.anahandles = opt.anahandles(i3(i2)); % seems like swapping the order
  opt.anahandles = opt.anahandles(:)';
  set(opt.anahandles, 'tag', 'ana');

  % for zooming purposes
  opt.axis = [opt.handlesaxes(1).XLim opt.handlesaxes(1).YLim opt.handlesaxes(1).ZLim];
  opt.redrawmarkers = true;
  opt.reinit = false; % do not redraw orthoplots

elseif opt.reinit
  ft_plot_ortho(opt.ana, 'transform', opt.mri{opt.currmri}.transform, 'coordsys', opt.mri{opt.currmri}.coordsys, 'location', opt.pos, 'style', 'subplot', 'surfhandle', opt.anahandles, 'update', opt.update, 'doscale', false, 'clim', opt.clim, 'unit', opt.mri{opt.currmri}.unit);

  fprintf('==================================================================================\n');
  lab = 'crosshair';
  switch opt.mri{opt.currmri}.unit
    case 'mm'
      fprintf('%10s at [%.1f %.1f %.1f] %s\n', lab, opt.pos, opt.mri{opt.currmri}.unit);
    case 'cm'
      fprintf('%10s at [%.2f %.2f %.2f] %s\n', lab, opt.pos, opt.mri{opt.currmri}.unit);
    case 'm'
      fprintf('%10s at [%.4f %.4f %.4f] %s\n', lab, opt.pos, opt.mri{opt.currmri}.unit);
    otherwise
      fprintf('%10s at [%f %f %f] %s\n', lab, opt.pos, opt.mri{opt.currmri}.unit);
  end
  opt.reinit = false; % do not redraw orthoplots
end

% zoom
xi = opt.pos(1);
yi = opt.pos(2);
zi = opt.pos(3);
xloadj = ((xi-opt.axis(1))-(xi-opt.axis(1))*opt.zoom);
xhiadj = ((opt.axis(2)-xi)-(opt.axis(2)-xi)*opt.zoom);
yloadj = ((yi-opt.axis(3))-(yi-opt.axis(3))*opt.zoom);
yhiadj = ((opt.axis(4)-yi)-(opt.axis(4)-yi)*opt.zoom);
zloadj = ((zi-opt.axis(5))-(zi-opt.axis(5))*opt.zoom);
zhiadj = ((opt.axis(6)-zi)-(opt.axis(6)-zi)*opt.zoom);
axis(h1, [xi-xloadj xi+xhiadj yi-yloadj yi+yhiadj zi-zloadj zi+zhiadj]);
axis(h2, [xi-xloadj xi+xhiadj yi-yloadj yi+yhiadj zi-zloadj zi+zhiadj]);
axis(h3, [xi-xloadj xi+xhiadj yi-yloadj yi+yhiadj]);

if opt.zoom>0
  % the coordsys labels fall outside the subplots when zoomed in
  delete(findall(h, 'Type', 'text', 'Tag', 'coordsyslabel_x'));
  delete(findall(h, 'Type', 'text', 'Tag', 'coordsyslabel_y'));
  delete(findall(h, 'Type', 'text', 'Tag', 'coordsyslabel_z'));
end

if opt.init
  % draw the crosshairs for the first time
  delete(findobj(opt.mainfig, 'Type', 'Line')); % get rid of old crosshairs (to facilitate switching scans)
  hch1 = ft_plot_crosshair([xi        yi-yloadj zi], 'parent', h1, 'color', 'yellow'); % was [xi 1 zi], now corrected for zoom
  hch2 = ft_plot_crosshair([xi+xhiadj yi        zi], 'parent', h2, 'color', 'yellow'); % was [opt.dim(1) yi zi], now corrected for zoom
  hch3 = ft_plot_crosshair([xi        yi        zi], 'parent', h3, 'color', 'yellow'); % was [xi yi opt.dim(3)], now corrected for zoom
  opt.handlescross  = [hch1(:)';hch2(:)';hch3(:)'];
else
  % update the existing crosshairs, don't change the handles
  ft_plot_crosshair([xi        yi-yloadj zi], 'handle', opt.handlescross(1, :));
  ft_plot_crosshair([xi+xhiadj yi        zi], 'handle', opt.handlescross(2, :));
  ft_plot_crosshair([xi        yi        zi], 'handle', opt.handlescross(3, :));
end

if opt.showcrosshair
  set(opt.handlescross, 'Visible', 'on');
else
  set(opt.handlescross, 'Visible', 'off');
end

% draw markers
if opt.showmarkers && opt.redrawmarkers
  delete(findobj(opt.mainfig, 'Tag', 'marker')); % remove previous markers
  delete(findobj(opt.mainfig, 'Tag', 'label')); % remove previous labels
  idx = find(~cellfun(@isempty,opt.markerlab)); % non-empty markers
  if ~isempty(idx)
    for i=1:numel(idx)
      markerlab_sel{i,1} = opt.markerlab{idx(i),1};
      markerpos_sel(i,:) = opt.markerpos{idx(i),1};
    end

    tmp1 = markerpos_sel(:,1);
    tmp2 = markerpos_sel(:,2);
    tmp3 = markerpos_sel(:,3);

    subplot(opt.handlesaxes(1));
    if ~opt.global % filter markers distant to the current slice (N units and further)
      posj_idx = find( abs(tmp2 - repmat(yi,size(tmp2))) < opt.markerdist);
      posi = tmp1(posj_idx);
      posj = tmp2(posj_idx);
      posk = tmp3(posj_idx);
    else % plot all markers on the current slice
      posj_idx = 1:numel(tmp1);
      posi = tmp1;
      posj = tmp2;
      posk = tmp3;
    end
    if ~isempty(posi)
      hold on
      plot3(posi, repmat(yi-yloadj,size(posj)), posk, 'Marker', '+', 'LineStyle', 'none', 'Color', 'r', 'Tag', 'marker'); % [xi yi-yloadj zi]
      if opt.showlabels
        for i=1:numel(posj_idx)
          text(posi(i), yi-yloadj, posk(i), markerlab_sel{posj_idx(i),1}, 'Color', [1 .5 0], 'Clipping', 'on', 'Tag', 'label');
        end
      end
      hold off
    end

    subplot(opt.handlesaxes(2));
    if ~opt.global % filter markers distant to the current slice (N units and further)
      posi_idx = find( abs(tmp1 - repmat(xi,size(tmp1))) < opt.markerdist);
      posi = tmp1(posi_idx);
      posj = tmp2(posi_idx);
      posk = tmp3(posi_idx);
    else % plot all markers on the current slice
      posi_idx = 1:numel(tmp1);
      posi = tmp1;
      posj = tmp2;
      posk = tmp3;
    end
    if ~isempty(posj)
      hold on
      plot3(repmat(xi+xhiadj,size(posi)), posj, posk, 'marker', '+', 'linestyle', 'none', 'color', 'r'); % [xi+xhiadj yi zi]
      if opt.showlabels
        for i=1:numel(posi_idx)
          text(posi(i)+xhiadj, posj(i), posk(i), markerlab_sel{posi_idx(i),1}, 'color', [1 .5 0], 'clipping', 'on');
        end
      end
      hold off
    end

    subplot(opt.handlesaxes(3));
    if ~opt.global % filter markers distant to the current slice (N units and further)
      posk_idx = find( abs(tmp3 - repmat(zi,size(tmp3))) < opt.markerdist);
      posi = tmp1(posk_idx);
      posj = tmp2(posk_idx);
      posk = tmp3(posk_idx);
    else % plot all markers on the current slice
      posk_idx = 1:numel(tmp1);
      posi = tmp1;
      posj = tmp2;
      posk = tmp3;
    end
    if ~isempty(posk)
      hold on
      plot3(posi, posj, repmat(zi,size(posk)), 'marker', '+', 'linestyle', 'none', 'color', 'r'); % [xi yi zi]
      if opt.showlabels
        for i=1:numel(posk_idx)
          text(posi(i), posj(i), zi, markerlab_sel{posk_idx(i),1}, 'color', [1 .5 0], 'clipping', 'on');
        end
      end
      hold off
    end
    clear markerlab_sel markerpos_sel
    opt.redrawmarkers = false;
    opt.redrawscattermarkers = true;
  end % if idx
end % if showmarkers

% also update the scatter appendix
if opt.scatter
  cb_redraw_scatter(h);
end

% make the last current axes current again
sel = findobj('type', 'axes', 'tag', tag);
if ~isempty(sel)
  set(opt.mainfig, 'currentaxes', sel(1));
end

% do not initialize on the next call
opt.init = false;
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw_scatter(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

if ~isfield(opt, 'scatterfig') % initiate in case the figure does not yet exist
  opt.scatterfig = figure(...
    'Name', [mfilename ' appendix'],...
    'Units', 'normalized', ...
    'Color', [1 1 1], ...
    'Visible', 'on');
  set(opt.scatterfig, 'CloseRequestFcn', @cb_scattercleanup);
  opt.scatterfig_h1 = axes('position', [0.02 0.02 0.96 0.96]);
  set(opt.scatterfig_h1, 'DataAspectRatio', get(opt.handlesaxes(1), 'DataAspectRatio'));

  dim = opt.mri{opt.currmri}.dim;
  transform = opt.mri{opt.currmri}.transform;
  % determine the corner points of the box that encompasses the whole data
  % extend the box with half a voxel in all directions to get the outer edge
  corner_vc = [
    0.5        0.5        0.5
    0.5+dim(1) 0.5        0.5
    0.5+dim(1) 0.5+dim(2) 0.5
    0.5        0.5+dim(2) 0.5
    0.5        0.5        0.5+dim(3)
    0.5+dim(1) 0.5        0.5+dim(3)
    0.5+dim(1) 0.5+dim(2) 0.5+dim(3)
    0.5        0.5+dim(2) 0.5+dim(3)
    ];
  corner_hc = ft_warp_apply(transform, corner_vc);
  % xlabel('x'); ylabel('y'); zlabel('z');
  axis equal; axis vis3d; axis off
  axis([min(corner_hc(:,1)) max(corner_hc(:,1)) min(corner_hc(:,2)) max(corner_hc(:,2)) min(corner_hc(:,3)) max(corner_hc(:,3))]);
  view(3);

  % scatter range sliders
  opt.scatterfig_h23text = uicontrol('Style', 'text',...
    'String', 'Intensity',...
    'Units', 'normalized', ...
    'Position', [.85+0.03 .26 .1 0.04],...
    'BackgroundColor', [1 1 1], ...
    'HandleVisibility', 'on');

  opt.scatterfig_h2 = uicontrol('Style', 'slider', ...
    'Parent', opt.scatterfig, ...
    'Min', 0, 'Max', 1, ...
    'Value', opt.slim(1), ...
    'Units', 'normalized', ...
    'Position', [.85+.02 .06 .05 .2], ...
    'Callback', @cb_scatterminslider);

  opt.scatterfig_h3 = uicontrol('Style', 'slider', ...
    'Parent', opt.scatterfig, ...
    'Min', 0, 'Max', 1, ...
    'Value', opt.slim(2), ...
    'Units', 'normalized', ...
    'Position', [.85+.07 .06 .05 .2], ...
    'Callback', @cb_scattermaxslider);

  hskullstrip = uicontrol('Style', 'togglebutton', ...
    'Parent', opt.scatterfig, ...
    'String', 'Skullstrip', ...
    'Value', 0, ...
    'Units', 'normalized', ...
    'Position', [.88 .88 .1 .1], ...
    'HandleVisibility', 'on', ...
    'Callback', @cb_skullstrip);

  % datacursor mode options
  opt.scatterfig_dcm = datacursormode(opt.scatterfig);
  set(opt.scatterfig_dcm, ...
    'DisplayStyle', 'datatip', ...
    'SnapToDataVertex', 'off', ...
    'Enable', 'on', ...
    'UpdateFcn', @cb_scatter_dcm);

  % draw the crosshair for the first time
  opt.handlescross2 = ft_plot_crosshair(opt.pos, 'parent', opt.scatterfig_h1, 'color', 'blue');

  % instructions to the user
  fprintf(strcat(...
    '5. Scatterplot viewing options:\n',...
    '   a. use the Data Cursor, Rotate 3D, Pan, and Zoom tools to navigate to electrodes in 3D space\n'));

  opt.redrawscatter = 1;
  opt.redrawscattermarkers = 1;
else
  set(0, 'CurrentFigure', opt.scatterfig) % make current figure (needed for ft_plot_crosshair)
end

if opt.redrawscatter
  delete(findobj(opt.scatterfig, 'Type', 'scatter')); % remove previous scatters
  msize = round(2000/opt.mri{opt.currmri}.dim(3)); % headsize (20 cm) / z slices
  inc = abs(opt.slim(2)-opt.slim(1))/4; % color increments
  for r = 1:4 % 4 color layers to encode peaks
    lim1 = opt.slim(1) + r*inc - inc;
    lim2 = opt.slim(1) + r*inc;
    voxind = find(opt.ana>lim1 & opt.ana<lim2);
    [x,y,z] = ind2sub(opt.mri{opt.currmri}.dim, voxind);
    pos = ft_warp_apply(opt.mri{opt.currmri}.transform, [x,y,z]);
    hold on 
    scatter3(opt.scatterfig_h1, pos(:,1), pos(:,2), pos(:,3), msize, 'Marker', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.8-(r*.2) .8-(r*.2) .8-(r*.2)], 'Tag', 'marker');
  end
  opt.redrawscatter = 0;
end

if opt.showmarkers && opt.redrawscattermarkers % plot the markers
  delete(findobj(opt.scatterfig, 'Type', 'line', 'Marker', '+')); % remove all scatterfig markers
  delete(findobj(opt.scatterfig, 'Type', 'text')); % remove all scatterfig labels
  idx = find(~cellfun(@isempty,opt.markerlab)); % non-empty markers
  if ~isempty(idx)
    for i=1:numel(idx)
      markerlab_sel{i,1} = opt.markerlab{idx(i),1};
      markerpos_sel(i,:) = opt.markerpos{idx(i),1};
    end
    plot3(opt.scatterfig_h1, markerpos_sel(:,1), markerpos_sel(:,2),markerpos_sel(:,3), 'marker', '+', 'linestyle', 'none', 'color', 'r'); % plot the markers
    if opt.showlabels
      for i=1:size(markerpos_sel,1)
        text(opt.scatterfig_h1, markerpos_sel(i,1), markerpos_sel(i,2), markerpos_sel(i,3), markerlab_sel{i,1}, 'color', [1 .5 0]);
      end
    end
    clear markerlab_sel markerpos_sel
  end % if idx
  opt.redrawscattermarkers = false;
end

% update the existing crosshairs, don't change the handles
ft_plot_crosshair(opt.pos, 'handle', opt.handlescross2);
if opt.showcrosshair
  set(opt.handlescross2, 'Visible', 'on');
else
  set(opt.handlescross2, 'Visible', 'off');
end

setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw_headshape(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

figure(h); % make current figure

delete(findobj(h, 'Type', 'line'));   % remove all lines 
delete(findobj(h, 'Type', 'Patch'));  % remove all markers
delete(findobj(h, 'Type', 'text'));   % remove all labels
delete(findall(h, 'Type', 'light'));  % remove all lights
delete(findobj(h, 'tag',  'headshape')); % remove the headshape

% plot the faces of the 2D or 3D triangulation
if ~opt.showsurface
  ft_plot_mesh(removefields(opt.headshape, 'color'), 'tag', 'headshape', 'facecolor', 'none', 'edgecolor', 'none', 'vertexcolor', 'none');
elseif isfield(opt.headshape, 'color') && opt.showcolors
  ft_plot_mesh(opt.headshape, 'tag', 'headshape');
else
  ft_plot_mesh(removefields(opt.headshape, 'color'), 'tag', 'headshape', 'facecolor', 'skin', 'edgecolor', 'none', 'facealpha', 1);
end

% apply some light and reflection for the 3D effect
lighting gouraud
if opt.camlight
  camlight
else
  uniformlight
end
if opt.shiny
  material shiny
else
  material dull
end
% alpha 0.9

if opt.showmarkers
  idx = find(~cellfun(@isempty,opt.markerlab)); % find the non-empty markers
  if ~isempty(idx)
    elec = keepfields(opt.headshape, {'unit', 'coordsys'});
    elec.elecpos = cat(1, opt.markerpos{idx});
    elec.label   = cat(1, opt.markerlab(idx));
    elec.elecori = elec.elecpos;
    elec.elecori(:,1) = elec.elecori(:,1) - mean(opt.headshape.pos(:,1));
    elec.elecori(:,2) = elec.elecori(:,2) - mean(opt.headshape.pos(:,2));
    elec.elecori(:,3) = elec.elecori(:,3) - mean(opt.headshape.pos(:,3));
    for i=1:numel(elec.label)
      elec.elecori(i,:) = elec.elecori(i,:) / norm(elec.elecori(i,:));
    end

    if opt.showlabels
      ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'label');
    else
      ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'off');
    end
  end % if not empty
end % if showmarkers

setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end

% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h, 'currentaxes');
tag     = get(curr_ax, 'tag');

if isempty(key)
  % this happens if you press the apple key
  key = '';
end

% the following code is largely shared by FT_SOURCEPLOT, FT_VOLUMEREALIGN, FT_INTERACTIVEREALIGN, FT_MESHREALIGN, FT_ELECTRODEPLACEMENT
if strcmp(opt.method, 'volume')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch key
    case {'' 'shift+shift' 'alt-alt' 'control+control' 'command-0'}
      % do nothing

    case 'h'
      cb_help(h)

    case 'q'
      setappdata(h, 'opt', opt);
      cb_quit(h);

    case '1'
      subplot(opt.handlesaxes(1));

    case '2'
      subplot(opt.handlesaxes(2));

    case '3'
      subplot(opt.handlesaxes(3));

    case 'g' % global/local elec view (h9) toggle
      if isequal(opt.global, 0)
        opt.global = 1;
        set(opt.handlesaxes(9), 'Value', 1);
      elseif isequal(opt.global, 1)
        opt.global = 0;
        set(opt.handlesaxes(9), 'Value', 0);
      end
      opt.redrawmarkers = true; % draw markers
      setappdata(h, 'opt', opt);
      cb_redraw(h);

    case 'l' % elec label view (h8) toggle
      if isequal(opt.showlabels, 0)
        opt.showlabels = 1;
        set(opt.handlesaxes(8), 'Value', 1);
      elseif isequal(opt.showlabels, 1)
        opt.showlabels = 0;
        set(opt.handlesaxes(8), 'Value', 0);
      end
      opt.redrawmarkers = true; % draw markers
      setappdata(h, 'opt', opt);
      cb_redraw(h);

    case 'm' % magnet (h7) toggle
      if isequal(opt.magnet, 0)
        opt.magnet = 1;
        set(opt.handlesaxes(7), 'Value', 1);
      elseif isequal(opt.magnet, 1)
        opt.magnet = 0;
        set(opt.handlesaxes(7), 'Value', 0);
      end
      setappdata(h, 'opt', opt);

    case {28 29 30 31 'leftarrow' 'rightarrow' 'uparrow' 'downarrow'}
      % update the view to a new position
      if     strcmp(tag, 'ik') && (strcmp(key, 'i') || strcmp(key, 'uparrow')    || isequal(key, 30)), opt.pos(3) = opt.pos(3)+1; opt.update = [1 1 1]; %[0 0 1];
      elseif strcmp(tag, 'ik') && (strcmp(key, 'j') || strcmp(key, 'leftarrow')  || isequal(key, 28)), opt.pos(1) = opt.pos(1)-1; opt.update = [1 1 1]; %[0 1 0];
      elseif strcmp(tag, 'ik') && (strcmp(key, 'k') || strcmp(key, 'rightarrow') || isequal(key, 29)), opt.pos(1) = opt.pos(1)+1; opt.update = [1 1 1]; %[0 1 0];
      elseif strcmp(tag, 'ik') && (strcmp(key, 'm') || strcmp(key, 'downarrow')  || isequal(key, 31)), opt.pos(3) = opt.pos(3)-1; opt.update = [1 1 1]; %[0 0 1];
      elseif strcmp(tag, 'ij') && (strcmp(key, 'i') || strcmp(key, 'uparrow')    || isequal(key, 30)), opt.pos(2) = opt.pos(2)+1; opt.update = [1 1 1]; %[1 0 0];
      elseif strcmp(tag, 'ij') && (strcmp(key, 'j') || strcmp(key, 'leftarrow')  || isequal(key, 28)), opt.pos(1) = opt.pos(1)-1; opt.update = [1 1 1]; %[0 1 0];
      elseif strcmp(tag, 'ij') && (strcmp(key, 'k') || strcmp(key, 'rightarrow') || isequal(key, 29)), opt.pos(1) = opt.pos(1)+1; opt.update = [1 1 1]; %[0 1 0];
      elseif strcmp(tag, 'ij') && (strcmp(key, 'm') || strcmp(key, 'downarrow')  || isequal(key, 31)), opt.pos(2) = opt.pos(2)-1; opt.update = [1 1 1]; %[1 0 0];
      elseif strcmp(tag, 'jk') && (strcmp(key, 'i') || strcmp(key, 'uparrow')    || isequal(key, 30)), opt.pos(3) = opt.pos(3)+1; opt.update = [1 1 1]; %[0 0 1];
      elseif strcmp(tag, 'jk') && (strcmp(key, 'j') || strcmp(key, 'leftarrow')  || isequal(key, 28)), opt.pos(2) = opt.pos(2)-1; opt.update = [1 1 1]; %[1 0 0];
      elseif strcmp(tag, 'jk') && (strcmp(key, 'k') || strcmp(key, 'rightarrow') || isequal(key, 29)), opt.pos(2) = opt.pos(2)+1; opt.update = [1 1 1]; %[1 0 0];
      elseif strcmp(tag, 'jk') && (strcmp(key, 'm') || strcmp(key, 'downarrow')  || isequal(key, 31)), opt.pos(3) = opt.pos(3)-1; opt.update = [1 1 1]; %[0 0 1];
      else
        % do nothing
      end
      opt.pos = min(opt.pos(:)', opt.axis([2 4 6])); % avoid out-of-bounds
      opt.pos = max(opt.pos(:)', opt.axis([1 3 5]));
      opt.reinit = true; % redraw orthoplots

      setappdata(h, 'opt', opt);
      cb_redraw(h);

      % contrast scaling
    case {43 'add' 'shift+equal'}  % + or numpad +
      if isempty(opt.clim)
        opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
      end
      % reduce color scale range by 10%
      cscalefactor = (opt.clim(2)-opt.clim(1))/10;
      %opt.clim(1) = opt.clim(1)+cscalefactor;
      opt.clim(2) = opt.clim(2)-cscalefactor;
      setappdata(h, 'opt', opt);
      cb_redraw(h);

    case {45 'subtract' 'hyphen' 'shift+hyphen'} % - or numpad -
      if isempty(opt.clim)
        opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
      end
      % increase color scale range by 10%
      cscalefactor = (opt.clim(2)-opt.clim(1))/10;
      %opt.clim(1) = opt.clim(1)-cscalefactor;
      opt.clim(2) = opt.clim(2)+cscalefactor;
      setappdata(h, 'opt', opt);
      cb_redraw(h);

    case 'c' % for crosshair
      opt.showcrosshair = ~opt.showcrosshair;
      setappdata(h, 'opt', opt);
      cb_redraw(h);

    case 'f' % for fiducials
      opt.showmarkers = ~opt.showmarkers;
      setappdata(h, 'opt', opt);
      cb_redraw(h);

    case 3 % right mouse click
      % add point to a list
      l1 = get(get(gca, 'xlabel'), 'string');
      l2 = get(get(gca, 'ylabel'), 'string');
      switch l1
        case 'i'
          xc = d1;
        case 'j'
          yc = d1;
        case 'k'
          zc = d1;
      end
      switch l2
        case 'i'
          xc = d2;
        case 'j'
          yc = d2;
        case 'k'
          zc = d2;
      end
      pnt = [pnt; xc yc zc];

    case 2 % middle mouse click
      l1 = get(get(gca, 'xlabel'), 'string');
      l2 = get(get(gca, 'ylabel'), 'string');

      % remove the previous point
      if size(pnt,1)>0
        pnt(end,:) = [];
      end

      if l1=='i' && l2=='j'
        updatepanel = [1 2 3];
      elseif l1=='i' && l2=='k'
        updatepanel = [2 3 1];
      elseif l1=='j' && l2=='k'
        updatepanel = [3 1 2];
      end

    case 'pageup'
      opt.zoom = opt.zoom + 0.1;
      opt.zoom = min(opt.zoom, 0.9); % between 0 and 0.9
      fprintf('increasing zoom to %f\n', opt.zoom);
      set(opt.handlesaxes(10), 'Value', opt.zoom);
      setappdata(h, 'opt', opt);
      cb_redraw(h);
      drawnow

    case 'pagedown'
      opt.zoom = opt.zoom - 0.1;
      opt.zoom = max(opt.zoom, 0); % between 0 and 0.9
      fprintf('decreasing zoom to %f\n', opt.zoom);
      set(opt.handlesaxes(10), 'Value', opt.zoom);
      setappdata(h, 'opt', opt);
      cb_redraw(h);
      drawnow

    otherwise
      % do nothing
  end % switch key

elseif strcmp(opt.method, 'headshape')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch key
    case {'' 'shift+shift' 'alt-alt' 'control+control' 'command-0'}
      % do nothing

    case 'h'
      cb_help(h)

    case 'q'
      setappdata(h, 'opt', opt);
      cb_quit(h);

    case 'v' % camlight angle reset
      delete(findall(h, 'Type', 'light')) % shut out the lights
      if opt.camlight
        camlight
      else
        uniformlight
      end

    otherwise
      % do nothing
  end % switch key
end % if method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonpress(h, eventdata)

h = getparent(h);
cb_getposition(h);
opt = getappdata(h, 'opt');
if strcmp(opt.method, 'volume') % only redraw volume/orthoplot
  switch get(h, 'selectiontype')
    case 'normal'
      % just update to new position, nothing else to be done here
      cb_redraw(h);
    case 'alt'
      set(h, 'windowbuttonmotionfcn', @cb_tracemouse);
      cb_redraw(h);
    otherwise
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonrelease(h, eventdata)

set(h, 'windowbuttonmotionfcn', '');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_tracemouse(h, eventdata)

h = getparent(h);
cb_getposition(h);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_getposition(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

if strcmp(opt.method, 'volume')
  curr_ax = get(h,       'currentaxes');
  tag = get(curr_ax, 'tag');
  if ~isempty(tag) && ~opt.init
    pos     = mean(get(curr_ax, 'currentpoint'));
    if strcmp(tag, 'ik')
      opt.pos([1 3])  = pos([1 3]);
      opt.update = [1 1 1];
    elseif strcmp(tag, 'ij')
      opt.pos([1 2])  = pos([1 2]);
      opt.update = [1 1 1];
    elseif strcmp(tag, 'jk')
      opt.pos([2 3])  = pos([2 3]);
      opt.update = [1 1 1];
    end
    opt.pos = min(opt.pos(:)', opt.axis([2 4 6])); % avoid out-of-bounds
    opt.pos = max(opt.pos(:)', opt.axis([1 3 5]));
  end
  if opt.magradius>0
    opt = magnetize(opt);
  end
  opt.reinit = true; % redraw orthoplots
  opt.redrawmarkers = true; % draw markers

elseif strcmp(opt.method, 'headshape')
  h2 = get(gca, 'children'); % get the object handles
  iscorrect = false(size(h2));
  for i=1:length(h2) % select the correct objects
    try
      pos = get(h2(i),'vertices');
      tri = get(h2(i),'faces');
      if ~isempty(opt.headshape) && isequal(opt.headshape.pos, pos) && isequal(opt.headshape.tri, tri)
        % it is the same object that the user has plotted before
        iscorrect(i) = true;
      elseif isempty(opt.headshape)
        % assume that it is the same object that the user has plotted before
        iscorrect(i) = true;
      end
    end
  end
  h2 = h2(iscorrect);
  opt.pos = select3d(h2)'; % enforce column direction
  if ~isempty(opt.pos)
    delete(findobj(h,'Type', 'Line', 'Marker', '+', 'Color', [0 0 0])) % remove previous crosshairs
    hold on; plot3(opt.pos(:,1), opt.pos(:,2), opt.pos(:,3), 'marker', '+', 'linestyle', 'none', 'color', [0 0 0], 'MarkerSize', 12); % plot the crosshair
  end
  %opt.pos = ft_select_point3d(opt.headshape, 'nearest', true, 'multiple', false, 'marker', '+'); % FIXME: this gets stuck in a loop waiting for any abritrary buttonpress
end
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quit(h, eventdata)

opt = getappdata(h, 'opt');
if isfield(opt, 'scatterfig')
  cb_scattercleanup(opt.scatterfig);
end
opt.quit = true;
setappdata(h, 'opt', opt);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_minslider(h4, eventdata)

newlim = get(h4, 'value');
h = getparent(h4);
opt = getappdata(h, 'opt');
opt.clim(1) = newlim;
opt.reinit = true; % redraw orthoplots
fprintf('==================================================================================\n');
fprintf(' contrast limits updated to [%.03f %.03f]\n', opt.clim);
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_maxslider(h5, eventdata)

newlim = get(h5, 'value');
h = getparent(h5);
opt = getappdata(h, 'opt');
opt.clim(2) = newlim;
opt.reinit = true; % redraw orthoplots
fprintf('==================================================================================\n');
fprintf(' contrast limits updated to [%.03f %.03f]\n', opt.clim);
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_intensityslider(h4, eventdata) % java intensity range slider - not fully functional

loval = get(h4, 'value');
hival = get(h4, 'highvalue');
h = getparent(h4); % this fails: The name 'parent' is not an accessible property for an instance of class 'com.jidesoft.swing.RangeSlider'.
opt = getappdata(h, 'opt');
opt.clim = [loval hival];
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_eleclistbox(h6, eventdata)

h = getparent(h6);
opt = getappdata(h, 'opt');

elecidx = get(h6, 'Value'); % chosen elec
listtopidx = get(h6, 'ListboxTop'); % ensure listbox does not move upon label selec
if ~isempty(elecidx)
  if numel(elecidx)>1
    fprintf('too many labels selected\n');
    return
  end
  if isempty(opt.pos)
    fprintf('no valid point was selected\n');
    return
  end
  eleclis = cellstr(get(h6, 'String')); % all labels
  eleclab = eleclis{elecidx}; % this elec's label

  % toggle electrode status and assign markers
  if contains(eleclab, 'silver') % not yet, check
    fprintf('==================================================================================\n');
    fprintf(' assigning marker %s\n', opt.label{elecidx,1});
    eleclab = regexprep(eleclab, '"silver"', '"black"'); % replace font color
    opt.markerlab{elecidx,1} = opt.label{elecidx,1}; % assign marker label
    opt.markerpos{elecidx,1} = opt.pos; % assign marker position
    opt.redrawmarkers = true; % draw markers
  elseif contains(eleclab, 'black') % already chosen before, move cursor to marker or uncheck
    if strcmp(get(h, 'SelectionType'), 'normal') % single click to move cursor to
      fprintf('==================================================================================\n');
      fprintf(' moving cursor to marker %s\n', opt.label{elecidx,1});
      opt.pos = opt.markerpos{elecidx,1}; % move cursor to marker position
      opt.reinit = true; % redraw orthoplots
    elseif strcmp(get(h, 'SelectionType'), 'open') % double click to uncheck
      fprintf('==================================================================================\n');
      fprintf(' removing marker %s\n', opt.label{elecidx,1});
      eleclab = regexprep(eleclab, '"black"', '"silver"'); % replace font color
      opt.markerlab{elecidx,1} = []; % assign marker label
      opt.markerpos{elecidx,1} = []; % assign marker position
      opt.redrawmarkers = true; % remove markers
    end
  end

  % update plot
  eleclis{elecidx} = eleclab;
  set(h6, 'String', eleclis);
  set(h6, 'ListboxTop', listtopidx); % ensure listbox does not move upon label selec
  setappdata(h, 'opt', opt);
  cb_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_magnetbutton(h7, eventdata)

h = getparent(h7);
opt = getappdata(h, 'opt');
radii = get(h7, 'String');
opt.magradius = str2double(radii{get(h7, 'value')});
fprintf('==================================================================================\n');
fprintf(' changed magnet radius to %.1f %s\n', opt.magradius, opt.mri{opt.currmri}.unit);
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = magnetize(opt)

try
  oldpos = opt.pos;
  vox = round(ft_warp_apply(inv(opt.mri{opt.currmri}.transform), oldpos)); % head to vox coord (for indexing within anatomy)
  xsel = vox(1)+(-opt.magradius:opt.magradius);
  ysel = vox(2)+(-opt.magradius:opt.magradius);
  zsel = vox(3)+(-opt.magradius:opt.magradius);
  cubic = opt.ana(xsel, ysel, zsel);
  if strcmp(opt.magtype, 'peak')
    % find the peak intensity voxel within the cube
    [val, idx] = max(cubic(:));
    [ix, iy, iz] = ind2sub(size(cubic), idx);
  elseif strcmp(opt.magtype, 'trough')
    % find the trough intensity voxel within the cube
    [val, idx] = min(cubic(:));
    [ix, iy, iz] = ind2sub(size(cubic), idx);
  elseif strcmp(opt.magtype, 'weighted')
    % find the weighted center of mass in the cube
    dim = size(cubic);
    [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
    cubic = cubic./sum(cubic(:));
    ix = (X(:)' * cubic(:));
    iy = (Y(:)' * cubic(:));
    iz = (Z(:)' * cubic(:));
  elseif strcmp(opt.magtype, 'peakweighted')
    % find the peak intensity voxel and then the center of mass
    [val, idx] = max(cubic(:));
    [ix, iy, iz] = ind2sub(size(cubic), idx);
    vox = [ix, iy, iz] + vox - opt.magradius - 1; % move cursor to peak
    xsel = vox(1)+(-opt.magradius:opt.magradius);
    ysel = vox(2)+(-opt.magradius:opt.magradius);
    zsel = vox(3)+(-opt.magradius:opt.magradius);
    cubic = opt.ana(xsel, ysel, zsel);
    dim = size(cubic);
    [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
    cubic = cubic./sum(cubic(:));
    ix = (X(:)' * cubic(:));
    iy = (Y(:)' * cubic(:));
    iz = (Z(:)' * cubic(:));
  elseif strcmp(opt.magtype, 'troughweighted')
    % find the peak intensity voxel and then the center of mass
    [val, idx] = min(cubic(:));
    [ix, iy, iz] = ind2sub(size(cubic), idx);
    vox = [ix, iy, iz] + vox - opt.magradius - 1; % move cursor to trough
    xsel = vox(1)+(-opt.magradius:opt.magradius);
    ysel = vox(2)+(-opt.magradius:opt.magradius);
    zsel = vox(3)+(-opt.magradius:opt.magradius);
    cubic = opt.ana(xsel, ysel, zsel);
    dim = size(cubic);
    [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
    cubic = 1-cubic;
    cubic = cubic./(sum(cubic(:)));
    ix = (X(:)' * cubic(:));
    iy = (Y(:)' * cubic(:));
    iz = (Z(:)' * cubic(:));
  end
  % adjust the indices for the selection
  voxadj = [ix, iy, iz] + vox - opt.magradius - 1;
  newpos = ft_warp_apply(opt.mri{opt.currmri}.transform, voxadj);
  if ~any(isnan(newpos))
    % this fails if the selection is at the edge of the volume
    fprintf('==================================================================================\n');
    fprintf(' clicked at [%.1f %.1f %.1f], %s magnetized adjustment [%.1f %.1f %.1f] %s\n', oldpos, opt.magtype, newpos-oldpos, opt.mri{opt.currmri}.unit);
    opt.pos = newpos;
  end
catch
  % this fails if the selection is at the edge of the volume
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_labelsbutton(h8, eventdata)

h = getparent(h8);
opt = getappdata(h, 'opt');
opt.showlabels = get(h8, 'value');
opt.redrawmarkers = true; % draw markers
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_colorsbutton(h9, eventdata)

h = getparent(h9);
opt = getappdata(h, 'opt');
opt.showcolors = get(h9, 'value');
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_camlight(h10, eventdata)

h = getparent(h10);
opt = getappdata(h, 'opt');
opt.camlight = get(h10, 'value');
setappdata(h, 'opt', opt);
% also apply it to the current figure
delete(findall(h, 'Type', 'light'));
if opt.camlight
  camlight
else
  uniformlight
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_shiny(h11, eventdata)

h = getparent(h11);
opt = getappdata(h, 'opt');
opt.shiny = get(h11, 'value');
setappdata(h, 'opt', opt);
% also apply it to the current figure
if opt.shiny
  material shiny
else
  material dull
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_surfacebutton(h9, eventdata)

h = getparent(h9);
opt = getappdata(h, 'opt');
opt.showsurface = get(h9, 'value');
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_globalbutton(h9, eventdata)

h = getparent(h9);
opt = getappdata(h, 'opt');
opt.global = get(h9, 'value');
opt.redrawmarkers = true; % draw markers
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_zoomslider(h10, eventdata)

h = getparent(h10);
opt = getappdata(h, 'opt');
opt.zoom = round(get(h10, 'value')*10)/10;
opt.redrawmarkers = true; % draw markers while ensuring they dont fall outside the bounding box
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scatterbutton(h11, eventdata)

h = getparent(h11);
opt = getappdata(h, 'opt');
opt.scatter = get(h11, 'value'); % update value
setappdata(h, 'opt', opt);
if isfield(opt, 'scatterfig') && ~opt.scatter % if already open but shouldn't, close it
  cb_scattercleanup(opt.scatterfig);
end
if opt.scatter
  cb_redraw_scatter(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scattercleanup(hObject, eventdata)

h = findobj('type', 'figure', 'name', mfilename);
opt = getappdata(h, 'opt');
opt.scatter = 0;
set(opt.handlesaxes(11), 'Value', 0);
opt = rmfield(opt, 'scatterfig');
setappdata(h, 'opt', opt);
delete(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scatterminslider(h2, eventdata)

h = findobj('type', 'figure', 'name', mfilename);
opt = getappdata(h, 'opt');
opt.slim(1) = get(h2, 'value');
fprintf('==================================================================================\n');
fprintf(' scatter limits updated to [%.03f %.03f]\n', opt.slim);
opt.redrawscatter = 1;
setappdata(h, 'opt', opt);
cb_redraw_scatter(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scattermaxslider(h3, eventdata)

h = findobj('type', 'figure', 'name', mfilename);
opt = getappdata(h, 'opt');
opt.slim(2) = get(h3, 'value');
fprintf('==================================================================================\n');
fprintf(' scatter limits updated to [%.03f %.03f]\n', opt.slim);
opt.redrawscatter = 1;
setappdata(h, 'opt', opt);
cb_redraw_scatter(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dcm_txt = cb_scatter_dcm(hObject, eventdata)

h = findobj('type', 'figure', 'name', mfilename);
opt = getappdata(h, 'opt');
opt.pos = get(eventdata, 'Position'); % current datamarker position
opt.reinit = true; % redraw orthoplots
dcm_txt = ''; % ['index = [' num2str(opt.pos) ']'];

if strcmp(get(opt.scatterfig_dcm, 'Enable'), 'on') % update appl and figures
  setappdata(h, 'opt', opt);
  cb_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_skullstrip(hObject, eventdata)

h = findobj('type', 'figure', 'name',mfilename);
opt = getappdata(h, 'opt');
if get(hObject, 'value') && ~isfield(opt.mri{opt.currmri}, 'dat_strip') % skullstrip
  fprintf('==================================================================================\n');
  fprintf(' stripping the skull - this may take a few minutes\n')
  tmp = keepfields(opt.mri{opt.currmri}, {'anatomy', 'dim', 'transform', 'unit', 'coordsys'});
  cfg = [];
  cfg.output = 'skullstrip';
  seg = ft_volumesegment(cfg, tmp);
  dmin = min(seg.anatomy(:));
  dmax = max(seg.anatomy(:));
  opt.mri{opt.currmri}.dat_strip = (seg.anatomy-dmin)./(dmax-dmin); % range between 0 and 1
  opt.ana = opt.mri{opt.currmri}.dat_strip; % overwrite with skullstrip
  clear seg tmp dmin dmax
elseif ~get(hObject, 'value') && isfield(opt.mri{opt.currmri}, 'dat_strip') % use original again
  opt.ana = opt.mri{opt.currmri}.dat;
elseif get(hObject, 'value') && isfield(opt.mri{opt.currmri}, 'dat_strip') % use skullstrip again
  opt.ana = opt.mri{opt.currmri}.dat_strip;
end
opt.redrawscatter = 1;
setappdata(h, 'opt', opt);
figure(h);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scanbutton(h12, eventdata)

h = getparent(h12);
opt = getappdata(h, 'opt');
opt.scan = get(h12, 'value');

opt.mri{opt.currmri}.clim = opt.clim; % store current scan's clim
opt.mri{opt.currmri}.slim = opt.slim; % store current scan's scatter clim

opt.currmri = opt.currmri+1; % rotate to next scan
if opt.currmri>numel(opt.mri)
  opt.currmri = 1; % restart at 1
end
opt.ana = opt.mri{opt.currmri}.dat;
opt.clim = opt.mri{opt.currmri}.clim; % use next scan's clim
set(opt.handlesaxes(4), 'Value', opt.clim(1)); % update minslider
set(opt.handlesaxes(5), 'Value', opt.clim(2)); % update maxslider
opt.slim = opt.mri{opt.currmri}.slim; % use next scan's scatter clim
opt.handlesaxes(1:3) = opt.mri{opt.currmri}.axes;
opt.init = true;

setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uniformlight
% place 6 lights uniformly around the object
l = lightangle(0,  90); set(l, 'Color', 0.45*[1 1 1])
l = lightangle(0, -90); set(l, 'Color', 0.45*[1 1 1])
l = lightangle(  0, 0); set(l, 'Color', 0.45*[1 1 1])
l = lightangle( 90, 0); set(l, 'Color', 0.45*[1 1 1])
l = lightangle(180, 0); set(l, 'Color', 0.45*[1 1 1])
l = lightangle(270, 0); set(l, 'Color', 0.45*[1 1 1])
