function [elec] = ft_electrodeplacement(cfg, varargin)

% FT_ELECTRODEPLACEMENT allows placing electrodes on a volume or headshape.
% The different methods are described in detail below.
%
% VOLUME - Navigate an orthographic display of a volume (e.g. CT or
% MR scan), and assign an electrode label to the current crosshair location
% by clicking on a label in the eletrode list. You can undo the selection by
% clicking on the same label again. The electrode labels shown in the list
% can be prespecified using cfg.channel when calling ft_electrodeplacement.
% The zoom slider allows zooming in at the location of the crosshair.
% The intensity sliders allow thresholding the image's low and high values.
% The magnet feature transports the crosshair to the nearest peak intensity
% voxel, within a certain voxel radius of the selected location.
% The labels feature displays the labels of the selected electrodes within
% the orthoplot. The global feature allows toggling the view between all
% and near-crosshair markers.
%
% HEADSHAPE - Navigate a triangulated head/brain surface, and assign
% an electrode location by clicking on the brain. The electrode
% is placed on the triangulation itself. FIXME: this needs updating
%
% Use as
%   [elec] = ft_electrodeplacement(cfg, mri)
% where the input mri should be an anatomical CT or MRI volume
% Use as
%   [elec] = ft_electrodeplacement(cfg, headshape)
% where the input headshape should be a surface triangulation
%
% The configuration can contain the following options
%   cfg.method         = string representing the method for aligning or placing the electrodes
%                        'mri'             place electrodes in a brain volume
%                        'headshape'       place electrodes on the head surface
%   cfg.parameter      = string, field in data (default = 'anatomy' if present in data)
%   cfg.channel        = Nx1 cell-array with selection of channels (default = '1','2', ...)
%   cfg.elec           = struct containing previously placed electrodes (this overwrites cfg.channel)
%   cfg.clim           = color range of the data (default = [0 1], i.e. the full range)
%   cfg.magtype        = string representing the 'magnet' type used for placing the electrodes
%                        'peak'            place electrodes at peak intensity voxel (default)
%                        'trough'          place electrodes at trough intensity voxel
%                        'weighted'        place electrodes at center-of-mass
%   cfg.magradius      = number representing the radius for the cfg.magtype based search (default = 2)
%
% See also FT_ELECTRODEREALIGN, FT_VOLUMEREALIGN

% Copyright (C) 2015, Arjen Stolk & Robert Oostenveld
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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% ensure that old and unsupported options are not being relied on by the end-user's script
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2837
cfg = ft_checkconfig(cfg, 'renamed', {'viewdim', 'axisratio'});

% set the defaults
cfg.method        = ft_getopt(cfg, 'method');                 % volume, headshape
cfg.parameter     = ft_getopt(cfg, 'parameter',    'anatomy');
cfg.channel       = ft_getopt(cfg, 'channel',             []); % default will be determined further down {'1', '2', ...}
cfg.elec          = ft_getopt(cfg, 'elec',                []); % use previously placed electrodes
cfg.renderer      = ft_getopt(cfg, 'renderer',      'opengl');
% view options
cfg.clim          = ft_getopt(cfg, 'clim',             [0 1]); % initial volume intensity limit voxels
cfg.markerdist    = ft_getopt(cfg, 'markerdist',           5); % marker-slice distance for including in the view
% magnet options
cfg.magtype       = ft_getopt(cfg, 'magtype',         'peak'); % detect peaks or troughs or center-of-mass
cfg.magradius     = ft_getopt(cfg, 'magradius',            2); % specify the physical unit radius
cfg.voxelratio    = ft_getopt(cfg, 'voxelratio',      'data'); % display size of the voxel, 'data' or 'square'
cfg.axisratio     = ft_getopt(cfg, 'axisratio',       'data'); % size of the axes of the three orthoplots, 'square', 'voxel', or 'data'


if isempty(cfg.method) && ~isempty(varargin)
  % the default determines on the input data
  switch ft_datatype(varargin{1})
    case 'volume'
      cfg.method = 'volume';
    case 'mesh'
      cfg.method = 'headshape';
  end
end

% check if the input data is valid for this function
switch cfg.method
  case 'volume'
    mri = ft_checkdata(varargin{1}, 'datatype', 'volume', 'feedback', 'yes');
  case  {'headshape'}
    headshape = fixpos(varargin{1});
    headshape = ft_determine_coordsys(headshape);
end

switch cfg.method
  case 'headshape'
    % give the user instructions
    disp('Use the mouse to click on the desired electrode positions');
    disp('Afterwards you may have to update the electrode labels');
    disp('Press "r" to delete the last point add');
    disp('Press "+/-" to zoom in/out');
    disp('Press "w/a/s/d" to rotate');
    disp('Press "q" when you are done');
    % open a figure
    figure;
    % plot the faces of the 2D or 3D triangulation
    
    if isfield(headshape, 'color');
      skin = 'none';
      ft_plot_mesh(headshape);
    else
      skin = [255 213 119]/255;
      ft_plot_mesh(headshape,'facecolor', skin,'EdgeColor','none','facealpha',1);
      lighting gouraud
      material shiny
      camlight
    end
    
    % rotate3d on
    xyz = ft_select_point3d(headshape, 'nearest', false, 'multiple', true, 'marker', '*');
    numelec = size(xyz, 1);

    % construct the output electrode structure
    elec = keepfields(headshape, {'unit', 'coordsys'});
    elec.elecpos = xyz;
    for i=1:numelec
      try
        elec.label{i} = cfg.channel{i,1};
      catch
        elec.label{i} = sprintf('%d', i);
      end
    end

  case 'volume'
    
    % start building the figure
    h = figure(...
      'MenuBar', 'none',...
      'Name', mfilename,...
      'Units', 'normalized', ...
      'Color', [1 1 1], ...
      'Visible', 'on');
    set(h, 'windowbuttondownfcn', @cb_buttonpress);
    set(h, 'windowbuttonupfcn',   @cb_buttonrelease);
    set(h, 'windowkeypressfcn',   @cb_keyboard);
    set(h, 'CloseRequestFcn',     @cb_cleanup);
    set(h, 'renderer', cfg.renderer);

    % axes settings
    if strcmp(cfg.axisratio, 'voxel')
      % determine the number of voxels to be plotted along each axis
      axlen1 = mri.dim(1);
      axlen2 = mri.dim(2);
      axlen3 = mri.dim(3);
    elseif strcmp(cfg.axisratio, 'data')
      % determine the length of the edges along each axis
      [cp_voxel, cp_head] = cornerpoints(mri.dim, mri.transform);
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
    h1size(1) = 0.82*axlen1/(axlen1 + axlen2);
    h1size(2) = 0.82*axlen3/(axlen2 + axlen3);
    h2size(1) = 0.82*axlen2/(axlen1 + axlen2);
    h2size(2) = 0.82*axlen3/(axlen2 + axlen3);
    h3size(1) = 0.82*axlen1/(axlen1 + axlen2);
    h3size(2) = 0.82*axlen2/(axlen2 + axlen3);

    if strcmp(cfg.voxelratio, 'square')
      voxlen1 = 1;
      voxlen2 = 1;
      voxlen3 = 1;
    elseif strcmp(cfg.voxelratio, 'data')
      % the size of the voxel is scaled with the data
      [cp_voxel, cp_head] = cornerpoints(mri.dim, mri.transform);
      voxlen1 = norm(cp_head(2,:)-cp_head(1,:))/norm(cp_voxel(2,:)-cp_voxel(1,:));
      voxlen2 = norm(cp_head(4,:)-cp_head(1,:))/norm(cp_voxel(4,:)-cp_voxel(1,:));
      voxlen3 = norm(cp_head(5,:)-cp_head(1,:))/norm(cp_voxel(5,:)-cp_voxel(1,:));
    end

    % axis handles will hold the anatomical functional if present, along with labels etc.
    h1 = axes('position',[0.06                0.06+0.06+h3size(2) h1size(1) h1size(2)]);
    h2 = axes('position',[0.06+0.06+h1size(1) 0.06+0.06+h3size(2) h2size(1) h2size(2)]);
    h3 = axes('position',[0.06                0.06                h3size(1) h3size(2)]);

    set(h1, 'Tag', 'ik', 'Visible', 'off', 'XAxisLocation', 'top');
    set(h2, 'Tag', 'jk', 'Visible', 'off', 'YAxisLocation', 'right'); % after rotating in ft_plot_ortho this becomes top
    set(h3, 'Tag', 'ij', 'Visible', 'off');

    set(h1, 'DataAspectRatio', 1./[voxlen1 voxlen2 voxlen3]);
    set(h2, 'DataAspectRatio', 1./[voxlen1 voxlen2 voxlen3]);
    set(h3, 'DataAspectRatio', 1./[voxlen1 voxlen2 voxlen3]);

    xc = round(mri.dim(1)/2); % start with center view
    yc = round(mri.dim(2)/2);
    zc = round(mri.dim(3)/2);

    dat = double(mri.(cfg.parameter));
    dmin = min(dat(:));
    dmax = max(dat(:));
    dat  = (dat-dmin)./(dmax-dmin); % range between 0 and 1

    % intensity range sliders
    h45text = uicontrol('Style', 'text',...
      'String','Intensity',...
      'Units', 'normalized', ...
      'Position',[2*h1size(1)+0.03 h3size(2)+0.03 h1size(1)/4 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on');

    h4 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 1, ...
      'Value', cfg.clim(1), ...
      'Units', 'normalized', ...
      'Position', [2*h1size(1)+0.02 0.15+h3size(2)/3 0.05 h3size(2)/2-0.05], ...
      'Callback', @cb_minslider);

    h5 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 1, ...
      'Value', cfg.clim(2), ...
      'Units', 'normalized', ...
      'Position', [2*h1size(1)+0.07 0.15+h3size(2)/3 0.05 h3size(2)/2-0.05], ...
      'Callback', @cb_maxslider);

    % java intensity range slider (dual-knob slider): the java component gives issues when wanting to
    % access the opt structure
    % [jRangeSlider] = com.jidesoft.swing.RangeSlider(0,1,cfg.clim(1),cfg.clim(2));  % min,max,low,high
    % [jRangeSlider, h4] = javacomponent(jRangeSlider, [], h);
    % set(h4, 'Units', 'normalized', 'Position', [0.05+h1size(1) 0.07 0.07 h3size(2)], 'Parent', h);
    % set(jRangeSlider, 'Orientation', 1, 'PaintTicks', true, 'PaintLabels', true, ...
    %     'Background', java.awt.Color.white, 'StateChangedCallback', @cb_intensityslider);

    % electrode listbox
    if ~isempty(cfg.elec) % re-use previously placed (cfg.elec) electrodes
      cfg.channel = []; % ensure cfg.channel is empty, for filling it up
      for e = 1:numel(cfg.elec.label)
        cfg.channel{e,1} = cfg.elec.label{e};
        chanstring{e} = ['<HTML><FONT color="black">' cfg.channel{e,1} '</FONT></HTML>']; % hmtl'ize

        markerlab{e,1} = cfg.elec.label{e};
        markerpos{e,1} = cfg.elec.elecpos(e,:);
      end
    else % otherwise use standard / prespecified (cfg.channel) electrode labels
      if isempty(cfg.channel)
        for c = 1:150
          cfg.channel{c,1} = sprintf('%d', c);
        end
      end
      for c = 1:numel(cfg.channel)
        chanstring{c} = ['<HTML><FONT color="silver">' cfg.channel{c,1} '</FONT></HTML>']; % hmtl'ize

        markerlab{c,1} = {};
        markerpos{c,1} = zeros(0,3);
      end
    end

    h6 = uicontrol('Style', 'listbox', ...
      'Parent', h, ...
      'Value', [], 'Min', 0, 'Max', numel(chanstring), ...
      'Units', 'normalized', ...
      'Position', [0.07+h1size(1)+0.05 0.07 h1size(1)/2 h3size(2)], ...
      'Callback', @cb_eleclistbox, ...
      'String', chanstring);

    % switches / radio buttons
    h7 = uicontrol('Style', 'radiobutton',...
      'Parent', h, ...
      'Value', 1, ...
      'String','Magnet',...
      'Units', 'normalized', ...
      'Position',[2*h1size(1) 0.22 h1size(1)/3 0.05],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_magnetbutton);

    h8 = uicontrol('Style', 'radiobutton',...
      'Parent', h, ...
      'Value', 0, ...
      'String','Labels',...
      'Units', 'normalized', ...
      'Position',[2*h1size(1) 0.17 h1size(1)/3 0.05],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_labelsbutton);

    h9 = uicontrol('Style', 'radiobutton',...
      'Parent', h, ...
      'Value', 0, ...
      'String','Global',...
      'Units', 'normalized', ...
      'Position',[2*h1size(1) 0.12 h1size(1)/3 0.05],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_globalbutton);

    hscatter = uicontrol('Style', 'radiobutton',...
      'Parent', h, ...
      'Value', 0, ...
      'String','Scatter',...
      'Units', 'normalized', ...
      'Position',[2*h1size(1) 0.07 h1size(1)/3 0.05],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_scatterbutton);
    
    % zoom slider
    h10text = uicontrol('Style', 'text',...
      'String','Zoom',...
      'Units', 'normalized', ...
      'Position',[1.8*h1size(1)+0.01 h3size(2)+0.03 h1size(1)/4 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on');

    h10 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 0.9, ...
      'Value', 0, ...
      'Units', 'normalized', ...
      'Position', [1.8*h1size(1)+0.02 0.15+h3size(2)/3 0.05 h3size(2)/2-0.05], ...
      'SliderStep', [.1 .1], ...
      'Callback', @cb_zoomslider);

    % instructions to the user
    fprintf(strcat(...
      '1. Viewing options:\n',...
      '   a. use the left mouse button to navigate the image, or\n',...
      '   b. use the arrow keys to increase or decrease the slice number by one\n',...
      '2. Placement options:\n',...
      '   a. click an electrode label in the list to assign the crosshair location, or\n',...
      '   b. doubleclick a previously assigned electrode label to remove its marker\n',...
      '3. To finalize, close the window or press q on the keyboard\n'));

    % create structure to be passed to gui
    opt               = [];
    opt.dim           = mri.dim;
    opt.ijk           = [xc yc zc];
    opt.h1size        = h1size;
    opt.h2size        = h2size;
    opt.h3size        = h3size;
    opt.handlesaxes   = [h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 hscatter];
    opt.handlesfigure = h;
    opt.handlesmarker = [];
    opt.quit          = false;
    opt.ana           = dat;
    opt.update        = [1 1 1];
    opt.init          = true;
    opt.tag           = 'ik';
    opt.mri           = mri;
    opt.showcrosshair = true;
    opt.vox           = [opt.ijk]; % voxel coordinates (physical units)
    opt.pos           = ft_warp_apply(mri.transform, opt.ijk); % head coordinates (e.g. mm)
    opt.showlabels    = 0;
    opt.label         = cfg.channel;
    opt.magnet        = get(h7, 'Value');
    opt.magradius     = cfg.magradius;
    opt.magtype       = cfg.magtype;
    opt.showmarkers   = true;
    opt.global        = get(h9, 'Value'); % show all markers in the current slices
    opt.scatter       = get(hscatter, 'Value'); % additional scatterplot
    opt.slim          = [.8 1]; % 80% - maximum
    opt.markerlab     = markerlab;
    opt.markerpos     = markerpos;
    opt.markerdist    = cfg.markerdist; % hidden option
    opt.clim          = cfg.clim;
    opt.zoom          = 0;
    if isfield(mri, 'unit') && ~strcmp(mri.unit, 'unknown')
      opt.unit = mri.unit;  % this is shown in the feedback on screen
    else
      opt.unit = '';        % this is not shown
    end

    setappdata(h, 'opt', opt);
    cb_redraw(h);

    while(opt.quit==0)
      uiwait(h);
      opt = getappdata(h, 'opt');
    end
    delete(h);

    % collect the results
    elec.label  = {};
    elec.elecpos = [];
    elec.chanpos = [];
    elec.tra = [];
    for i=1:length(opt.markerlab)
      if ~isempty(opt.markerlab{i,1})
        elec.label = [elec.label; opt.markerlab{i,1}];
        elec.elecpos = [elec.elecpos; opt.markerpos{i,1}];
      end
    end
    elec.chanpos  = elec.elecpos; % identicial to elecpos
    elec.tra = eye(size(elec.elecpos,1));
    if isfield(mri, 'unit')
      elec.unit  = mri.unit;
    end
    if isfield(mri, 'coordsys')
      elec.coordsys = mri.coordsys;
    end
    
end % switch method

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   mri
ft_postamble provenance elec
ft_postamble history    elec
ft_postamble savevar    elec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)
tic
h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h, 'currentaxes');
tag = get(curr_ax, 'tag');

mri = opt.mri;

h1 = opt.handlesaxes(1);
h2 = opt.handlesaxes(2);
h3 = opt.handlesaxes(3);

xi = opt.ijk(1);
yi = opt.ijk(2);
zi = opt.ijk(3);

if any([xi yi zi] > mri.dim) || any([xi yi zi] <= 0)
  return;
end

opt.ijk = [xi yi zi 1]';
opt.ijk = opt.ijk(1:3)';

% construct a string with user feedback
str1 = sprintf('voxel %d, index [%d %d %d]', sub2ind(mri.dim(1:3), xi, yi, zi), opt.ijk);

if opt.init
  ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', opt.ijk, 'style', 'subplot', 'parents', [h1 h2 h3], 'update', opt.update, 'doscale', false,'clim', opt.clim);

  opt.anahandles = findobj(opt.handlesfigure, 'type', 'surface')';
  parenttag = get(opt.anahandles,'parent');
  parenttag{1} = get(parenttag{1}, 'tag');
  parenttag{2} = get(parenttag{2}, 'tag');
  parenttag{3} = get(parenttag{3}, 'tag');
  [i1,i2,i3] = intersect(parenttag, {'ik';'jk';'ij'});
  opt.anahandles = opt.anahandles(i3(i2)); % seems like swapping the order
  opt.anahandles = opt.anahandles(:)';
  set(opt.anahandles, 'tag', 'ana');

  % for zooming purposes
  opt.axis = zeros(1,6);
  opt.axis([1 3 5]) = 0.5;
  opt.axis([2 4 6]) = size(opt.ana) + 0.5;
else
  ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', opt.ijk, 'style', 'subplot', 'surfhandle', opt.anahandles, 'update', opt.update, 'doscale', false,'clim', opt.clim);

  if all(round([xi yi zi])<=mri.dim) && all(round([xi yi zi])>0)
    fprintf('==================================================================================\n');
    str = sprintf('voxel %d, index [%d %d %d]', sub2ind(mri.dim(1:3), round(xi), round(yi), round(zi)), round([xi yi zi]));

    lab = 'crosshair';
    opt.vox = [xi yi zi];
    ind = sub2ind(mri.dim(1:3), round(opt.vox(1)), round(opt.vox(2)), round(opt.vox(3)));
    opt.pos = ft_warp_apply(mri.transform, opt.vox);
    switch opt.unit
      case 'mm'
        fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.1f %.1f %.1f] %s\n', lab, ind, opt.vox, opt.pos, opt.unit);
      case 'cm'
        fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.2f %.2f %.2f] %s\n', lab, ind, opt.vox, opt.pos, opt.unit);
      case 'm'
        fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.4f %.4f %.4f] %s\n', lab, ind, opt.vox, opt.pos, opt.unit);
      otherwise
        fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%f %f %f] %s\n', lab, ind, opt.vox, opt.pos, opt.unit);
    end
  end

end

% make the last current axes current again
sel = findobj('type','axes','tag',tag);
if ~isempty(sel)
  set(opt.handlesfigure, 'currentaxes', sel(1));
end

% zoom
xloadj = round((xi-opt.axis(1))-(xi-opt.axis(1))*opt.zoom);
xhiadj = round((opt.axis(2)-xi)-(opt.axis(2)-xi)*opt.zoom);
yloadj = round((yi-opt.axis(3))-(yi-opt.axis(3))*opt.zoom);
yhiadj = round((opt.axis(4)-yi)-(opt.axis(4)-yi)*opt.zoom);
zloadj = round((zi-opt.axis(5))-(zi-opt.axis(5))*opt.zoom);
zhiadj = round((opt.axis(6)-zi)-(opt.axis(6)-zi)*opt.zoom);
axis(h1, [xi-xloadj xi+xhiadj yi-yloadj yi+yhiadj zi-zloadj zi+zhiadj]);
axis(h2, [xi-xloadj xi+xhiadj yi-yloadj yi+yhiadj zi-zloadj zi+zhiadj]);
axis(h3, [xi-xloadj xi+xhiadj yi-yloadj yi+yhiadj]);

if opt.init
  % draw the crosshairs for the first time
  hch1 = crosshair([xi yi-yloadj zi], 'parent', h1, 'color', 'yellow'); % was [xi 1 zi], now corrected for zoom
  hch2 = crosshair([xi+xhiadj yi zi], 'parent', h2, 'color', 'yellow'); % was [opt.dim(1) yi zi], now corrected for zoom
  hch3 = crosshair([xi yi zi], 'parent', h3, 'color', 'yellow'); % was [xi yi opt.dim(3)], now corrected for zoom
  opt.handlescross  = [hch1(:)';hch2(:)';hch3(:)'];
  opt.handlesmarker = [];
else
  % update the existing crosshairs, don't change the handles
  crosshair([xi yi-yloadj zi], 'handle', opt.handlescross(1, :));
  crosshair([xi+xhiadj yi zi], 'handle', opt.handlescross(2, :));
  crosshair([xi yi zi], 'handle', opt.handlescross(3, :));
end

if opt.showcrosshair
  set(opt.handlescross,'Visible','on');
else
  set(opt.handlescross,'Visible','off');
end

delete(opt.handlesmarker(opt.handlesmarker(:)>0));
opt.handlesmarker = [];

% draw markers
idx = find(~cellfun(@isempty,opt.markerlab)); % non-empty markers
if ~isempty(idx)
  for i=1:numel(idx)
    markerlab{i,1} = opt.markerlab{idx(i),1};
    markerpos(i,:) = opt.markerpos{idx(i),1};
  end

  opt.vox2 = round(ft_warp_apply(inv(mri.transform), markerpos)); % head to vox
  tmp1 = opt.vox2(:,1);
  tmp2 = opt.vox2(:,2);
  tmp3 = opt.vox2(:,3);

  subplot(h1);
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
    opt.handlesmarker(:,1) = plot3(posi, repmat(yi-yloadj,size(posj)), posk, 'marker', '+', 'linestyle', 'none', 'color', 'r'); % [xi yi-yloadj zi]
    if opt.showlabels
      for i=1:numel(posj_idx)
        opt.handlesmarker(i,4) = text(posi(i), yi-yloadj, posk(i), markerlab{posj_idx(i),1}, 'color', 'b');
      end
    end
    hold off
  end

  subplot(h2);
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
    opt.handlesmarker(:,2) = plot3(repmat(xi+xhiadj,size(posi)), posj, posk, 'marker', '+', 'linestyle', 'none', 'color', 'r'); % [xi+xhiadj yi zi]
    if opt.showlabels
      for i=1:numel(posi_idx)
        opt.handlesmarker(i,5) = text(posi(i)+xhiadj, posj(i), posk(i), markerlab{posi_idx(i),1}, 'color', 'b');
      end
    end
    hold off
  end

  subplot(h3);
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
    opt.handlesmarker(:,3) = plot3(posi, posj, repmat(zi,size(posk)), 'marker', '+', 'linestyle', 'none', 'color', 'r'); % [xi yi zi]
    if opt.showlabels
      for i=1:numel(posk_idx)
        opt.handlesmarker(i,6) = text(posi(i), posj(i), zi, markerlab{posk_idx(i),1}, 'color', 'b');
      end
    end
    hold off
  end
end % for all markers

if isfield(opt, 'scatterfig')
  cb_scatterredraw(h); % also update the appendix
  figure(h); % FIXME: ugly as it switches forth and back to mainfig
end

% do not initialize on the next call
opt.init = false;
setappdata(h, 'opt', opt);
set(h, 'currentaxes', curr_ax);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scatterredraw(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

if opt.scatter % radiobutton on
  if ~isfield(opt, 'scatterfig') % if the figure does not yet exist, initiate
    opt.scatterfig = figure(...
      'Name', [mfilename ' appendix'],...
      'Units', 'normalized', ...
      'Color', [1 1 1], ...
      'Visible', 'on');
    set(opt.scatterfig, 'CloseRequestFcn', @cb_scattercleanup);
    opt.scatterfig_h1 = axes('position',[0.06 0.06 0.74 0.88]);
    set(opt.scatterfig_h1, 'DataAspectRatio', get(opt.handlesaxes(1), 'DataAspectRatio'));
    axis image;
    xlabel('x'); ylabel('y'); zlabel('z');
    
    % scatter range sliders
    opt.scatterfig_h23text = uicontrol('Style', 'text',...
      'String','Treshold',...
      'Units', 'normalized', ...
      'Position',[.85+0.03 .26 .1 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on');
    
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
    
    msize = round(2500/opt.mri.dim(3)); % headsize (25 cm) / z slices
    inc = abs(opt.slim(2)-opt.slim(1))/4; % color increments
    for r = 1:4 % 4 color layers to encode peaks
      lim1 = opt.slim(1) + r*inc - inc;
      lim2 = opt.slim(1) + r*inc;
      voxind = find(opt.ana>lim1 & opt.ana<lim2);
      [x,y,z] = ind2sub(opt.mri.dim, voxind);
      hold on; scatter3(x,y,z,msize,'Marker','s','MarkerEdgeColor','none','MarkerFaceColor',[.8-(r*.2) .8-(r*.2) .8-(r*.2)]);
    end
    
    % draw the crosshair for the first time
    opt.handlescross2 = crosshair([opt.ijk], 'parent', opt.scatterfig_h1, 'color', 'blue');
  end
  
  figure(opt.scatterfig);
  
  % update the existing crosshairs, don't change the handles
  crosshair([opt.ijk], 'handle', opt.handlescross2);
  if opt.showcrosshair
    set(opt.handlescross,'Visible','on');
  else
    set(opt.handlescross,'Visible','off');
  end
  
  % plot the markers
  if isfield(opt, 'vox2')
    delete(findobj(opt.scatterfig,'Type','line','Marker','+')); % remove previous markers
    plot3(opt.vox2(:,1),opt.vox2(:,2),opt.vox2(:,3), 'marker', '+', 'linestyle', 'none', 'color', 'r');
  end
end
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
  key = parseKeyboardEvent(eventdata);
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

% the following code is largely shared with FT_SOURCEPLOT
switch key
  case {'' 'shift+shift' 'alt-alt' 'control+control' 'command-0'}
    % do nothing

  case '1'
    subplot(opt.handlesaxes(1));

  case '2'
    subplot(opt.handlesaxes(2));

  case '3'
    subplot(opt.handlesaxes(3));

  case 'q'
    setappdata(h, 'opt', opt);
    cb_cleanup(h);

  case 'g' % global/local elec view (h9) toggle
    if isequal(opt.global, 0)
      opt.global = 1;
      set(opt.handlesaxes(9), 'Value', 1);
    elseif isequal(opt.global, 1)
      opt.global = 0;
      set(opt.handlesaxes(9), 'Value', 0);
    end
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
    if     strcmp(tag,'ik') && (strcmp(key,'i') || strcmp(key,'uparrow')    || isequal(key, 30)), opt.ijk(3) = opt.ijk(3)+1; opt.update = [0 0 1];
    elseif strcmp(tag,'ik') && (strcmp(key,'j') || strcmp(key,'leftarrow')  || isequal(key, 28)), opt.ijk(1) = opt.ijk(1)-1; opt.update = [0 1 0];
    elseif strcmp(tag,'ik') && (strcmp(key,'k') || strcmp(key,'rightarrow') || isequal(key, 29)), opt.ijk(1) = opt.ijk(1)+1; opt.update = [0 1 0];
    elseif strcmp(tag,'ik') && (strcmp(key,'m') || strcmp(key,'downarrow')  || isequal(key, 31)), opt.ijk(3) = opt.ijk(3)-1; opt.update = [0 0 1];
    elseif strcmp(tag,'ij') && (strcmp(key,'i') || strcmp(key,'uparrow')    || isequal(key, 30)), opt.ijk(2) = opt.ijk(2)+1; opt.update = [1 0 0];
    elseif strcmp(tag,'ij') && (strcmp(key,'j') || strcmp(key,'leftarrow')  || isequal(key, 28)), opt.ijk(1) = opt.ijk(1)-1; opt.update = [0 1 0];
    elseif strcmp(tag,'ij') && (strcmp(key,'k') || strcmp(key,'rightarrow') || isequal(key, 29)), opt.ijk(1) = opt.ijk(1)+1; opt.update = [0 1 0];
    elseif strcmp(tag,'ij') && (strcmp(key,'m') || strcmp(key,'downarrow')  || isequal(key, 31)), opt.ijk(2) = opt.ijk(2)-1; opt.update = [1 0 0];
    elseif strcmp(tag,'jk') && (strcmp(key,'i') || strcmp(key,'uparrow')    || isequal(key, 30)), opt.ijk(3) = opt.ijk(3)+1; opt.update = [0 0 1];
    elseif strcmp(tag,'jk') && (strcmp(key,'j') || strcmp(key,'leftarrow')  || isequal(key, 28)), opt.ijk(2) = opt.ijk(2)-1; opt.update = [1 0 0];
    elseif strcmp(tag,'jk') && (strcmp(key,'k') || strcmp(key,'rightarrow') || isequal(key, 29)), opt.ijk(2) = opt.ijk(2)+1; opt.update = [1 0 0];
    elseif strcmp(tag,'jk') && (strcmp(key,'m') || strcmp(key,'downarrow')  || isequal(key, 31)), opt.ijk(3) = opt.ijk(3)-1; opt.update = [0 0 1];
    else
      % do nothing
    end;

    setappdata(h, 'opt', opt);
    cb_redraw(h);

    % contrast scaling
  case {43 'shift+equal'}  % numpad +
    if isempty(opt.clim)
      opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
    end
    % reduce color scale range by 10%
    cscalefactor = (opt.clim(2)-opt.clim(1))/10;
    %opt.clim(1) = opt.clim(1)+cscalefactor;
    opt.clim(2) = opt.clim(2)-cscalefactor;
    setappdata(h, 'opt', opt);
    cb_redraw(h);

  case {45 'shift+hyphen'} % numpad -
    if isempty(opt.clim)
      opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
    end
    % increase color scale range by 10%
    cscalefactor = (opt.clim(2)-opt.clim(1))/10;
    %opt.clim(1) = opt.clim(1)-cscalefactor;
    opt.clim(2) = opt.clim(2)+cscalefactor;
    setappdata(h, 'opt', opt);
    cb_redraw(h);

  case 99  % 'c'
    opt.showcrosshair = ~opt.showcrosshair;
    setappdata(h, 'opt', opt);
    cb_redraw(h);

  case 102 % 'f'
    opt.showmarkers = ~opt.showmarkers;
    setappdata(h, 'opt', opt);
    cb_redraw(h);

  case 3 % right mouse click
    % add point to a list
    l1 = get(get(gca, 'xlabel'), 'string');
    l2 = get(get(gca, 'ylabel'), 'string');
    switch l1,
      case 'i'
        xc = d1;
      case 'j'
        yc = d1;
      case 'k'
        zc = d1;
    end
    switch l2,
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

  otherwise
    % do nothing

end % switch key

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonpress(h, eventdata)

h = getparent(h);
cb_getposition(h);
switch get(h, 'selectiontype')
  case 'normal'
    % just update to new position, nothing else to be done here
    cb_redraw(h);
  case 'alt'
    set(h, 'windowbuttonmotionfcn', @cb_tracemouse);
    cb_redraw(h);
  otherwise
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
curr_ax = get(h,       'currentaxes');
pos     = mean(get(curr_ax, 'currentpoint'));
tag = get(curr_ax, 'tag');
if ~isempty(tag) && ~opt.init
  if strcmp(tag, 'ik')
    opt.ijk([1 3])  = round(pos([1 3]));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'ij')
    opt.ijk([1 2])  = round(pos([1 2]));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'jk')
    opt.ijk([2 3])  = round(pos([2 3]));
    opt.update = [1 1 1];
  end
end
opt.ijk = min(opt.ijk(:)', opt.dim);
opt.ijk = max(opt.ijk(:)', [1 1 1]);

if opt.magnet % magnetize
  try
    center = opt.ijk;
    radius = opt.magradius;
    % FIXME here it would be possible to adjust the selection at the edges of the volume
    xsel = center(1)+(-radius:radius);
    ysel = center(2)+(-radius:radius);
    zsel = center(3)+(-radius:radius);
    cubic  = opt.ana(xsel, ysel, zsel);
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
      ix = round(X(:)' * cubic(:));
      iy = round(Y(:)' * cubic(:));
      iz = round(Z(:)' * cubic(:));
    end
    % adjust the indices for the selection
    opt.ijk = [ix, iy, iz] + center - radius - 1;
    fprintf('==================================================================================\n');
    fprintf(' clicked at [%d %d %d], %s magnetized adjustment [%d %d %d]\n', center, opt.magtype, opt.ijk-center);
  catch
    % this fails if the selection is at the edge of the volume
    warning('cannot magnetize at the edge of the volume');
  end
end % if magnetize
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_cleanup(h, eventdata)

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
function key = parseKeyboardEvent(eventdata)

key = eventdata.Key;
% handle possible numpad events (different for Windows and UNIX systems)
% NOTE: shift+numpad number does not work on UNIX, since the shift
% modifier is always sent for numpad events
if isunix()
  shiftInd = match_str(eventdata.Modifier, 'shift');
  if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
    % now we now it was a numpad keystroke (numeric character sent AND
    % shift modifier present)
    key = eventdata.Character;
    eventdata.Modifier(shiftInd) = []; % strip the shift modifier
  end
elseif ispc()
  if strfind(eventdata.Key, 'numpad')
    key = eventdata.Character;d
  end
end

if ~isempty(eventdata.Modifier)
  key = [eventdata.Modifier{1} '+' key];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_minslider(h4, eventdata)

newlim = get(h4, 'value');
h = getparent(h4);
opt = getappdata(h, 'opt');
opt.clim(1) = newlim;
fprintf('contrast limits updated to [%.03f %.03f]\n', opt.clim);
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
fprintf('contrast limits updated to [%.03f %.03f]\n', opt.clim);
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

elecidx = get(h6, 'Value'); % chosen elec
listtopidx = get(h6, 'ListboxTop'); % ensure listbox does not move upon label selec
if ~isempty(elecidx)
  if numel(elecidx)>1
    fprintf('too many labels selected\n');
    return
  end
  eleclis = cellstr(get(h6, 'String')); % all labels
  eleclab = eleclis{elecidx}; % this elec's label

  h = getparent(h6);
  opt = getappdata(h, 'opt');

  % toggle electrode status and assign markers
  if strfind(eleclab, 'silver') % not yet, check
    fprintf('assigning marker %s\n', opt.label{elecidx,1});
    eleclab = regexprep(eleclab, '"silver"','"black"'); % replace font color
    opt.markerlab{elecidx,1} = opt.label(elecidx,1); % assign marker label
    opt.markerpos{elecidx,1} = opt.pos; % assign marker position
  elseif strfind(eleclab, 'black') % already chosen before, move cusor to marker or uncheck
    if strcmp(get(h,'SelectionType'),'normal') % single click to move cursor to
      fprintf('moving cursor to marker %s\n', opt.label{elecidx,1});
      opt.ijk = ft_warp_apply(inv(opt.mri.transform), opt.markerpos{elecidx,1}); % move cursor to marker position
    elseif strcmp(get(h,'SelectionType'),'open') % double click to uncheck
      fprintf('removing marker %s\n', opt.label{elecidx,1});
      eleclab = regexprep(eleclab, '"black"','"silver"'); % replace font color
      opt.markerlab{elecidx,1} = {}; % assign marker label
      opt.markerpos{elecidx,1} = zeros(0,3); % assign marker position
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
opt.magnet = get(h7, 'value');
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_labelsbutton(h8, eventdata)

h = getparent(h8);
opt = getappdata(h, 'opt');
opt.showlabels = get(h8, 'value');
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_globalbutton(h9, eventdata)

h = getparent(h9);
opt = getappdata(h, 'opt');
opt.global = get(h9, 'value');
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_zoomslider(h10, eventdata)

h = getparent(h10);
opt = getappdata(h, 'opt');
opt.zoom = round(get(h10, 'value')*10)/10;
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scatterbutton(hscatter, eventdata)

h = getparent(hscatter);
opt = getappdata(h, 'opt');
opt.scatter = get(hscatter, 'value'); % update value
setappdata(h, 'opt', opt);
if isfield(opt, 'scatterfig') && ~opt.scatter % if already open but shouldn't, close it
  cb_scattercleanup(opt.scatterfig);
end
if opt.scatter
  cb_scatterredraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scattercleanup(hObject, eventdata)

h = findobj('type','figure','name',mfilename);
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

h = findobj('type','figure','name',mfilename);
opt = getappdata(h, 'opt');
opt.slim(1) = get(h2, 'value');
fprintf('scatter limits updated to [%.03f %.03f]\n', opt.slim);
setappdata(h, 'opt', opt);

delete(findobj('type','scatter')); % remove previous scatters
msize = round(2500/opt.mri.dim(3)); % headsize (25 cm) / z slices
inc = abs(opt.slim(2)-opt.slim(1))/4; % color increments
for r = 1:4 % 4 color layers to encode peaks
  lim1 = opt.slim(1) + r*inc - inc;
  lim2 = opt.slim(1) + r*inc;
  voxind = find(opt.ana>lim1 & opt.ana<lim2);
  [x,y,z] = ind2sub(opt.mri.dim, voxind);
  hold on; scatter3(x,y,z,msize,'Marker','s','MarkerEdgeColor','none','MarkerFaceColor',[.8-(r*.2) .8-(r*.2) .8-(r*.2)]);
end
cb_scatterredraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scattermaxslider(h3, eventdata)

h = findobj('type','figure','name',mfilename);
opt = getappdata(h, 'opt');
opt.slim(2) = get(h3, 'value');
fprintf('scatter limits updated to [%.03f %.03f]\n', opt.slim);
setappdata(h, 'opt', opt);

delete(findobj('type','scatter')); % remove previous scatters
msize = round(2500/opt.mri.dim(3)); % headsize (25 cm) / z slices
inc = abs(opt.slim(2)-opt.slim(1))/4; % color increments
for r = 1:4 % 4 color layers to encode peaks
  lim1 = opt.slim(1) + r*inc - inc;
  lim2 = opt.slim(1) + r*inc;
  voxind = find(opt.ana>lim1 & opt.ana<lim2);
  [x,y,z] = ind2sub(opt.mri.dim, voxind);
  hold on; scatter3(x,y,z,msize,'Marker','s','MarkerEdgeColor','none','MarkerFaceColor',[.8-(r*.2) .8-(r*.2) .8-(r*.2)]);
end
cb_scatterredraw(h);