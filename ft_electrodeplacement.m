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
% Use as
%   [elec] = ft_electrodeplacement(cfg, ct)
%   [elec] = ft_electrodeplacement(cfg, ct, mri)
% where the input mri should be an anatomical CT or MRI volume, or
%   [elec] = ft_electrodeplacement(cfg, headshape)
% where the input headshape should be a surface triangulation.
%
% The configuration can contain the following options
%   cfg.method         = string representing the method for placing the electrodes
%                        'mri'             interactively locate electrodes in a MRI or CT scan
%                        'headshape'       interactively locate electrodes on a head surface
%                        '1020'            automatically place electrodes on a head surface
%
% The following options apply to the mri method
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
% The following options apply to the 1020 method
%   cfg.fiducial.nas   = 1x3 vector with coordinates
%   cfg.fiducial.ini   = 1x3 vector with coordinates
%   cfg.fiducial.lpa   = 1x3 vector with coordinates
%   cfg.fiducial.rpa   = 1x3 vector with coordinates
%   cfg.feedback       = string, can be 'yes' or 'no' for detailled feedback (default = 'yes')
%
% See also FT_ELECTRODEREALIGN, FT_VOLUMEREALIGN, FT_VOLUMESEGMENT, FT_PREPARE_MESH

% Copyright (C) 2015-2017, Arjen Stolk, Sandon Griffin & Robert Oostenveld
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
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'mri', 'volume'});

% set the defaults
cfg.method        = ft_getopt(cfg, 'method');                  % volume, headshape, 1020
cfg.feedback      = ft_getopt(cfg, 'feedback',         'yes');
cfg.parameter     = ft_getopt(cfg, 'parameter',    'anatomy');
cfg.channel       = ft_getopt(cfg, 'channel',             []); % default will be determined further down {'1', '2', ...}
cfg.elec          = ft_getopt(cfg, 'elec',                []); % use previously placed electrodes
cfg.renderer      = ft_getopt(cfg, 'renderer',      'opengl');
% view options
cfg.clim          = ft_getopt(cfg, 'clim',             [0 1]); % initial volume intensity limit voxels
cfg.markerdist    = ft_getopt(cfg, 'markerdist',           5); % marker-slice distance view when ~global
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
    for v = 1:numel(varargin)
      mri{v} = ft_checkdata(varargin{v}, 'datatype', 'volume', 'feedback', 'yes', 'hascoordsys', 'yes', 'hasunit', 'yes');
    end
  case  {'headshape', '1020'}
    headshape = fixpos(varargin{1});
    headshape = ft_checkdata(headshape, 'hascoordsys', 'yes');
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
      h1 = axes('position',[0.02                0.02+0.04+h3size(2) h1size(1) h1size(2)]); % x z
      h2 = axes('position',[0.02+0.04+h1size(1) 0.02+0.04+h3size(2) h2size(1) h2size(2)]); % y z
      h3 = axes('position',[0.02                0.02                h3size(1) h3size(2)]); % x y
      
      set(h1, 'Tag', 'ik', 'Visible', 'off', 'XAxisLocation', 'top'); axis(h1,'equal');
      set(h2, 'Tag', 'jk', 'Visible', 'off', 'YAxisLocation', 'right'); axis(h2,'equal'); % after rotating in ft_plot_ortho this becomes top
      set(h3, 'Tag', 'ij', 'Visible', 'off'); axis(h3,'equal');
      
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
      'String','Intensity',...
      'Units', 'normalized', ...
      'Position',[2*mri{1}.h1size(1)-0.02 mri{1}.h3size(2)-0.02 mri{1}.h1size(1)/4 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on');
    
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
    chanlabel = {}; chanstring = {};
    markerlab = {}; markerpos = {};
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
          
          markerlab{end+1,1} = {};
          markerpos{end+1,1} = zeros(0,3);
        end
      end
    end
    if isempty(cfg.elec) && isempty(cfg.channel) % create electrode labels on-the-fly
      for c = 1:150
        chanlabel{end+1,1} = sprintf('%d', c);
        chanstring{end+1} = ['<HTML><FONT color="silver">' sprintf('%d', c) '</FONT></HTML>']; % hmtl'ize
        
        markerlab{end+1,1} = {};
        markerpos{end+1,1} = zeros(0,3);
      end
    end
    
    h6 = uicontrol('Style', 'listbox', ...
      'Parent', h, ...
      'Value', [], 'Min', 0, 'Max', numel(chanstring), ...
      'Units', 'normalized', ...
      'FontSize', 12, ...
      'Position', [mri{1}.h1size(1)+0.07 0.02 mri{1}.h2size(1)/2.5 mri{1}.h3size(2)], ...
      'Callback', @cb_eleclistbox, ...
      'String', chanstring);
        
    % switches / radio buttons
    h7 = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 1, ...
      'String','Magnet',...
      'Units', 'normalized', ...
      'Position',[2*mri{1}.h1size(1)-0.05 0.18 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_magnetbutton);
    
    h8 = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 0, ...
      'String','Labels',...
      'Units', 'normalized', ...
      'Position',[2*mri{1}.h1size(1)-0.05 0.14 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_labelsbutton);
    
    h9 = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 0, ...
      'String','Global',...
      'Units', 'normalized', ...
      'Position',[2*mri{1}.h1size(1)-0.05 0.10 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_globalbutton);
    
    hscatter = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 0, ...
      'String','Scatter',...
      'Units', 'normalized', ...
      'Position',[2*mri{1}.h1size(1)-0.05 0.06 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_scatterbutton);
    
    hscan = uicontrol('Style', 'checkbox',...
      'Parent', h, ...
      'Value', 0, ...
      'String','CT/MRI',...
      'Units', 'normalized', ...
      'Position',[2*mri{1}.h1size(1)-0.05 0.02 mri{1}.h1size(1)/3 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Visible', 'off', ...
      'Callback', @cb_scanbutton);
    if numel(mri)>1; set(hscan, 'Visible', 'on'); end % only when two scans are given as input
    
    % zoom slider
    h10text = uicontrol('Style', 'text',...
      'String','Zoom',...
      'Units', 'normalized', ...
      'Position',[1.8*mri{1}.h1size(1)-0.04 mri{1}.h3size(2)-0.02 mri{1}.h1size(1)/4 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on');
    
    h10 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 0.9, ...
      'Value', 0, ...
      'Units', 'normalized', ...
      'Position', [1.8*mri{1}.h1size(1)-0.03 0.10+mri{1}.h3size(2)/3 0.05 mri{1}.h3size(2)/2-0.05], ...
      'SliderStep', [.1 .1], ...
      'Callback', @cb_zoomslider);
    
    % instructions to the user
    fprintf(strcat(...
      '1. Orthoplot viewing options:\n',...
      '   a. use the left mouse button to navigate the image, or\n',...
      '   b. use the arrow keys to increase or decrease the slice number by one\n',...
      '2. Orthoplot placement options:\n',...
      '   a. click an electrode label in the list to assign it to the crosshair location, or\n',...
      '   b. doubleclick a previously assigned electrode label to remove its marker\n',...
      '3. To finalize, close the window or press q on the keyboard\n', ...
      '4. See Stolk, Griffin et al. (2017) for further electrode processing options\n'));
    
    % create structure to be passed to gui
    opt               = [];
    opt.label         = chanlabel;
    opt.axes          = [mri{1}.axes(1) mri{1}.axes(2) mri{1}.axes(3) h4 h5 h6 h7 h8 h9 h10 hscatter hscan];
    opt.mainfig       = h;
    opt.quit          = false;
    opt.update        = [1 1 1];
    opt.init          = true;
    opt.tag           = 'ik';
    opt.ana           = mri{1}.dat; % the plotted anatomy
    opt.mri           = mri;
    opt.showcrosshair = true;
    opt.pos           = [0 0 0]; % middle of the scan, head coordinates
    opt.showlabels    = false;
    opt.magnet        = get(h7, 'Value');
    opt.magradius     = cfg.magradius;
    opt.magtype       = cfg.magtype;
    opt.showmarkers   = true;
    opt.global        = get(h9, 'Value'); % show all markers in the current slices
    opt.scatter       = get(hscatter, 'Value'); % additional scatterplot
    opt.scan          = get(hscan, 'Value'); % switch scans
    opt.slim          = [.9 1]; % 90% - maximum
    opt.markerlab     = markerlab;
    opt.markerpos     = markerpos;
    opt.markerdist    = cfg.markerdist; % hidden option
    opt.clim          = cfg.clim;
    opt.zoom          = 0;
    
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
    while(opt.quit==0)
      uiwait(h);
      opt = getappdata(h, 'opt');
    end
    delete(h);
    
    % collect the results
    elec = keepfields(mri{1}, {'unit', 'coordsys'});
    elec.label   = {};
    elec.elecpos = [];
    for i=1:length(opt.markerlab)
      if ~isempty(opt.markerlab{i,1})
        elec.label = [elec.label; opt.markerlab{i,1}];
        elec.elecpos = [elec.elecpos; opt.markerpos{i,1}];
      end
    end
    elec.chanpos = elec.elecpos;
    elec.tra = eye(size(elec.elecpos,1));
    
  case '1020'
    
    % the placement procedure fails if the fiducials coincide with vertices
    dist = @(x, y) sqrt(sum(bsxfun(@minus, x, y).^2,2));
    tolerance = 0.1 * ft_scalingfactor('mm', headshape.unit);  % 0.1 mm
    nas = cfg.fiducial.nas;
    ini = cfg.fiducial.ini;
    lpa = cfg.fiducial.lpa;
    rpa = cfg.fiducial.rpa;
    if any(dist(headshape.pos, nas)<tolerance)
      warning('Nasion coincides with headshape vertex, addding random displacement of about %f %s', tolerance, headshape.unit);
      nas = nas + tolerance*randn(1,3);
    end
    if any(dist(headshape.pos, ini)<tolerance)
      warning('Inion coincides with headshape vertex, addding random displacement of about %f %s', tolerance, headshape.unit);
      ini = ini + tolerance*randn(1,3);
    end
    if any(dist(headshape.pos, lpa)<tolerance)
      warning('LPA coincides with headshape vertex, addding random displacement of about %f %s', tolerance, headshape.unit);
      lpa = lpa + tolerance*randn(1,3);
    end
    if any(dist(headshape.pos, rpa)<tolerance)
      warning('RPA coincides with headshape vertex, addding random displacement of about %f %s', tolerance, headshape.unit);
      rpa = rpa + tolerance*randn(1,3);
    end
    
    % place the electrodes automatically according to the fiducials
    [pos, lab] = elec1020_locate(headshape.pos, headshape.tri, nas, ini, lpa, rpa, istrue(cfg.feedback));
    % construct the output
    elec = keepfields(headshape, {'unit', 'coordsys'});
    elec.elecpos = pos;
    elec.label   = lab(:);
    
  otherwise
    error('unsupported method ''%s''', cfg.method);
    
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

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h, 'currentaxes');
tag = get(curr_ax, 'tag');

h1 = opt.axes(1);
h2 = opt.axes(2);
h3 = opt.axes(3);

if opt.init
  delete(findobj(opt.mainfig,'Type','Surface')); % get rid of old orthos (to facilitate switching scans)
  ft_plot_ortho(opt.ana, 'transform', opt.mri{1}.transform, 'location', opt.pos, 'style', 'subplot', 'parents', [h1 h2 h3], 'update', opt.update, 'doscale', false,'clim', opt.clim, 'unit', opt.mri{1}.unit);
  
  opt.anahandles = findobj(opt.mainfig,'Type','Surface')';
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
  opt.axis = [opt.axes(1).XLim opt.axes(1).YLim opt.axes(1).ZLim];
else
  ft_plot_ortho(opt.ana, 'transform', opt.mri{1}.transform, 'location', opt.pos, 'style', 'subplot', 'surfhandle', opt.anahandles, 'update', opt.update, 'doscale', false,'clim', opt.clim, 'unit', opt.mri{1}.unit);
  
  fprintf('==================================================================================\n');
  lab = 'crosshair';
  switch opt.mri{1}.unit
    case 'mm'
      fprintf('%10s at [%.1f %.1f %.1f] %s\n', lab, opt.pos, opt.mri{1}.unit);
    case 'cm'
      fprintf('%10s at [%.2f %.2f %.2f] %s\n', lab, opt.pos, opt.mri{1}.unit);
    case 'm'
      fprintf('%10s at [%.4f %.4f %.4f] %s\n', lab, opt.pos, opt.mri{1}.unit);
    otherwise
      fprintf('%10s at [%f %f %f] %s\n', lab, opt.pos, opt.mri{1}.unit);
  end
end

% make the last current axes current again
sel = findobj('type','axes','tag',tag);
if ~isempty(sel)
  set(opt.mainfig, 'currentaxes', sel(1));
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

if opt.init
  % draw the crosshairs for the first time
  delete(findobj(opt.mainfig,'Type','Line')); % get rid of old crosshairs (to facilitate switching scans)
  hch1 = ft_plot_crosshair([xi yi-yloadj zi], 'parent', h1, 'color', 'yellow'); % was [xi 1 zi], now corrected for zoom
  hch2 = ft_plot_crosshair([xi+xhiadj yi zi], 'parent', h2, 'color', 'yellow'); % was [opt.dim(1) yi zi], now corrected for zoom
  hch3 = ft_plot_crosshair([xi yi zi], 'parent', h3, 'color', 'yellow'); % was [xi yi opt.dim(3)], now corrected for zoom
  opt.handlescross  = [hch1(:)';hch2(:)';hch3(:)'];
  opt.redrawmarker = 1;
else
  % update the existing crosshairs, don't change the handles
  ft_plot_crosshair([xi yi-yloadj zi], 'handle', opt.handlescross(1, :));
  ft_plot_crosshair([xi+xhiadj yi zi], 'handle', opt.handlescross(2, :));
  ft_plot_crosshair([xi yi zi], 'handle', opt.handlescross(3, :));
end

if opt.showcrosshair
  set(opt.handlescross,'Visible','on');
else
  set(opt.handlescross,'Visible','off');
end

% draw markers
delete(findobj(h,'Type','Line','Marker','+')); % remove previous markers
delete(findobj(h,'Type','text')); % remove previous labels
if opt.showmarkers
  idx = find(~cellfun(@isempty,opt.markerlab)); % non-empty markers
  if ~isempty(idx)
    for i=1:numel(idx)
      opt.markerlab_sel{i,1} = opt.markerlab{idx(i),1};
      opt.markerpos_sel(i,:) = opt.markerpos{idx(i),1};
    end
    
    tmp1 = opt.markerpos_sel(:,1);
    tmp2 = opt.markerpos_sel(:,2);
    tmp3 = opt.markerpos_sel(:,3);
    
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
      plot3(posi, repmat(yi-yloadj,size(posj)), posk, 'marker', '+', 'linestyle', 'none', 'color', 'r'); % [xi yi-yloadj zi]
      if opt.showlabels
        for i=1:numel(posj_idx)
          text(posi(i), yi-yloadj, posk(i), opt.markerlab_sel{posj_idx(i),1}, 'color', [1 .5 0]);
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
      plot3(repmat(xi+xhiadj,size(posi)), posj, posk, 'marker', '+', 'linestyle', 'none', 'color', 'r'); % [xi+xhiadj yi zi]
      if opt.showlabels
        for i=1:numel(posi_idx)
          text(posi(i)+xhiadj, posj(i), posk(i), opt.markerlab_sel{posi_idx(i),1}, 'color', [1 .5 0]);
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
      plot3(posi, posj, repmat(zi,size(posk)), 'marker', '+', 'linestyle', 'none', 'color', 'r'); % [xi yi zi]
      if opt.showlabels
        for i=1:numel(posk_idx)
          text(posi(i), posj(i), zi, opt.markerlab_sel{posk_idx(i),1}, 'color', [1 .5 0]);
        end
      end
      hold off
    end
  end % for all markers
end % if showmarkers

% do not initialize on the next call
opt.init = false;
setappdata(h, 'opt', opt);
set(h, 'currentaxes', curr_ax);

% also update the appendix
if isfield(opt, 'scatterfig')
  cb_scatterredraw(h);
  figure(h); % this switches back to mainfig
end

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
    opt.scatterfig_h1 = axes('position',[0.02 0.02 0.96 0.96]);
    set(opt.scatterfig_h1, 'DataAspectRatio', get(opt.axes(1), 'DataAspectRatio'));
    axis square; axis tight; axis off; view([0 0]);
    xlabel('x'); ylabel('y'); zlabel('z');
    
    % scatter range sliders
    opt.scatterfig_h23text = uicontrol('Style', 'text',...
      'String','Intensity',...
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
      'Enable', 'off', ...
      'UpdateFcn', @cb_scatter_dcm);
    
    % draw the crosshair for the first time
    opt.handlescross2 = ft_plot_crosshair(opt.pos, 'parent', opt.scatterfig_h1, 'color', 'blue');
    
    % instructions to the user
    fprintf(strcat(...
      '5. Scatterplot viewing options:\n',...
      '   a. use the Data Cursor, Rotate 3D, Pan, and Zoom tools to navigate to electrodes in 3D space\n'));
    
    opt.redrawscatter = 1;
    opt.redrawmarker = 1;
  end
  
  figure(opt.scatterfig); % make current figure
  
  if opt.redrawscatter
    delete(findobj('type','scatter')); % remove previous scatters
    msize = round(2000/opt.mri{1}.dim(3)); % headsize (20 cm) / z slices
    inc = abs(opt.slim(2)-opt.slim(1))/4; % color increments
    for r = 1:4 % 4 color layers to encode peaks
      lim1 = opt.slim(1) + r*inc - inc;
      lim2 = opt.slim(1) + r*inc;
      voxind = find(opt.ana>lim1 & opt.ana<lim2);
      [x,y,z] = ind2sub(opt.mri{1}.dim, voxind);
      pos = ft_warp_apply(opt.mri{1}.transform, [x,y,z]);
      hold on; scatter3(pos(:,1),pos(:,2),pos(:,3),msize,'Marker','s','MarkerEdgeColor','none','MarkerFaceColor',[.8-(r*.2) .8-(r*.2) .8-(r*.2)]);
    end
    opt.redrawscatter = 0;
  end
  
  if opt.redrawmarker
    delete(findobj(opt.scatterfig,'Type','line','Marker','+')); % remove all scatterfig markers
    delete(findobj(opt.scatterfig,'Type','text')); % remove all scatterfig labels
    if opt.showmarkers && isfield(opt, 'markerpos_sel') % plot the markers
      plot3(opt.markerpos_sel(:,1),opt.markerpos_sel(:,2),opt.markerpos_sel(:,3), 'marker', '+', 'linestyle', 'none', 'color', 'r'); % plot the markers
      if opt.showlabels
        for i=1:size(opt.markerpos_sel,1)
          text(opt.markerpos_sel(i,1), opt.markerpos_sel(i,2), opt.markerpos_sel(i,3), opt.markerlab_sel{i,1}, 'color', [1 .5 0]);
        end
      end
    end
    opt.redrawmarker = 0;
  end

  % update the existing crosshairs, don't change the handles
  ft_plot_crosshair([opt.pos], 'handle', opt.handlescross2);
  if opt.showcrosshair
    set(opt.handlescross2,'Visible','on');
  else
    set(opt.handlescross2,'Visible','off');
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
    subplot(opt.axes(1));
    
  case '2'
    subplot(opt.axes(2));
    
  case '3'
    subplot(opt.axes(3));
    
  case 'q'
    setappdata(h, 'opt', opt);
    cb_cleanup(h);
    
  case 'g' % global/local elec view (h9) toggle
    if isequal(opt.global, 0)
      opt.global = 1;
      set(opt.axes(9), 'Value', 1);
    elseif isequal(opt.global, 1)
      opt.global = 0;
      set(opt.axes(9), 'Value', 0);
    end
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
  case 'l' % elec label view (h8) toggle
    if isequal(opt.showlabels, 0)
      opt.showlabels = 1;
      set(opt.axes(8), 'Value', 1);
    elseif isequal(opt.showlabels, 1)
      opt.showlabels = 0;
      set(opt.axes(8), 'Value', 0);
    end
    opt.redrawmarker = 1;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
  case 'm' % magnet (h7) toggle
    if isequal(opt.magnet, 0)
      opt.magnet = 1;
      set(opt.axes(7), 'Value', 1);
    elseif isequal(opt.magnet, 1)
      opt.magnet = 0;
      set(opt.axes(7), 'Value', 0);
    end
    setappdata(h, 'opt', opt);
    
  case {28 29 30 31 'leftarrow' 'rightarrow' 'uparrow' 'downarrow'}
    % update the view to a new position
    if     strcmp(tag,'ik') && (strcmp(key,'i') || strcmp(key,'uparrow')    || isequal(key, 30)), opt.pos(3) = opt.pos(3)+1; opt.update = [0 0 1];
    elseif strcmp(tag,'ik') && (strcmp(key,'j') || strcmp(key,'leftarrow')  || isequal(key, 28)), opt.pos(1) = opt.pos(1)-1; opt.update = [0 1 0];
    elseif strcmp(tag,'ik') && (strcmp(key,'k') || strcmp(key,'rightarrow') || isequal(key, 29)), opt.pos(1) = opt.pos(1)+1; opt.update = [0 1 0];
    elseif strcmp(tag,'ik') && (strcmp(key,'m') || strcmp(key,'downarrow')  || isequal(key, 31)), opt.pos(3) = opt.pos(3)-1; opt.update = [0 0 1];
    elseif strcmp(tag,'ij') && (strcmp(key,'i') || strcmp(key,'uparrow')    || isequal(key, 30)), opt.pos(2) = opt.pos(2)+1; opt.update = [1 0 0];
    elseif strcmp(tag,'ij') && (strcmp(key,'j') || strcmp(key,'leftarrow')  || isequal(key, 28)), opt.pos(1) = opt.pos(1)-1; opt.update = [0 1 0];
    elseif strcmp(tag,'ij') && (strcmp(key,'k') || strcmp(key,'rightarrow') || isequal(key, 29)), opt.pos(1) = opt.pos(1)+1; opt.update = [0 1 0];
    elseif strcmp(tag,'ij') && (strcmp(key,'m') || strcmp(key,'downarrow')  || isequal(key, 31)), opt.pos(2) = opt.pos(2)-1; opt.update = [1 0 0];
    elseif strcmp(tag,'jk') && (strcmp(key,'i') || strcmp(key,'uparrow')    || isequal(key, 30)), opt.pos(3) = opt.pos(3)+1; opt.update = [0 0 1];
    elseif strcmp(tag,'jk') && (strcmp(key,'j') || strcmp(key,'leftarrow')  || isequal(key, 28)), opt.pos(2) = opt.pos(2)-1; opt.update = [1 0 0];
    elseif strcmp(tag,'jk') && (strcmp(key,'k') || strcmp(key,'rightarrow') || isequal(key, 29)), opt.pos(2) = opt.pos(2)+1; opt.update = [1 0 0];
    elseif strcmp(tag,'jk') && (strcmp(key,'m') || strcmp(key,'downarrow')  || isequal(key, 31)), opt.pos(3) = opt.pos(3)-1; opt.update = [0 0 1];
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
    
  case 102 % 'f' for fiducials
    opt.showmarkers = ~opt.showmarkers;
    opt.redrawmarker = 1;
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

if opt.magnet % magnetize
  opt = magnetize(opt);
end
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
      opt.pos = opt.markerpos{elecidx,1}; % move cursor to marker position
    elseif strcmp(get(h,'SelectionType'),'open') % double click to uncheck
      fprintf('removing marker %s\n', opt.label{elecidx,1});
      eleclab = regexprep(eleclab, '"black"','"silver"'); % replace font color
      opt.markerlab{elecidx,1} = {}; % assign marker label
      opt.markerpos{elecidx,1} = []; % assign marker position
    end
  end
  
  % update plot
  eleclis{elecidx} = eleclab;
  set(h6, 'String', eleclis);
  set(h6, 'ListboxTop', listtopidx); % ensure listbox does not move upon label selec
  opt.redrawmarker = 1;
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
function opt = magnetize(opt)

try
  pos = opt.pos;
  vox = round(ft_warp_apply(inv(opt.mri{1}.transform), pos)); % head to vox coord (for indexing within anatomy)
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
    ix = round(X(:)' * cubic(:));
    iy = round(Y(:)' * cubic(:));
    iz = round(Z(:)' * cubic(:));
  end
  % adjust the indices for the selection
  voxadj = [ix, iy, iz] + vox - opt.magradius - 1;
  opt.pos = ft_warp_apply(opt.mri{1}.transform, voxadj);
  fprintf('==================================================================================\n');
  fprintf(' clicked at [%.1f %.1f %.1f], %s magnetized adjustment [%.1f %.1f %.1f] %s\n', pos, opt.magtype, opt.pos-pos, opt.mri{1}.unit);
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
opt.redrawmarker = 1;
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
set(opt.axes(11), 'Value', 0);
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
opt.redrawscatter = 1;
setappdata(h, 'opt', opt);
cb_scatterredraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scattermaxslider(h3, eventdata)

h = findobj('type','figure','name',mfilename);
opt = getappdata(h, 'opt');
opt.slim(2) = get(h3, 'value');
fprintf('scatter limits updated to [%.03f %.03f]\n', opt.slim);
opt.redrawscatter = 1;
setappdata(h, 'opt', opt);
cb_scatterredraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dcm_txt = cb_scatter_dcm(hObject, eventdata)

h = findobj('type','figure','name',mfilename);
opt = getappdata(h, 'opt');
opt.pos = get(eventdata, 'Position'); % current datamarker position
if opt.magnet % magnetize
  opt = magnetize(opt);
end
dcm_txt = ['']; % ['index = [' num2str(opt.pos) ']'];

if strcmp(get(opt.scatterfig_dcm, 'Enable'), 'on') % update appl and figures
  setappdata(h, 'opt', opt);
  figure(h);
  cb_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_skullstrip(hObject, eventdata)

h = findobj('type','figure','name',mfilename);
opt = getappdata(h, 'opt');
if get(hObject, 'value') && ~isfield(opt.mri{1}, 'dat_strip') % skullstrip
  fprintf('stripping the skull - this could take a few minutes\n')
  tmp = keepfields(opt.mri{1}, {'anatomy','dim','coordsys','unit','transform'});
  cfg = [];
  cfg.output = 'skullstrip';
  seg = ft_volumesegment(cfg, tmp);
  dmin = min(seg.anatomy(:));
  dmax = max(seg.anatomy(:));
  opt.mri{1}.dat_strip = (seg.anatomy-dmin)./(dmax-dmin); % range between 0 and 1 
  opt.ana = opt.mri{1}.dat_strip; % overwrite with skullstrip
  clear seg tmp dmin dmax
elseif ~get(hObject, 'value') && isfield(opt.mri{1}, 'dat_strip') % use original again
  opt.ana = opt.mri{1}.dat;
elseif get(hObject, 'value') && isfield(opt.mri{1}, 'dat_strip') % use skullstrip again
  opt.ana = opt.mri{1}.dat_strip;
end
opt.redrawscatter = 1;
setappdata(h, 'opt', opt);
figure(h);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_scanbutton(hscan, eventdata)

h = getparent(hscan);
opt = getappdata(h, 'opt');
opt.scan = get(hscan, 'value');

opt.mri{1}.clim = opt.clim; % store current clim
opt.mri{1}.slim = opt.slim; % store current scatter clim
opt.mri = flip(opt.mri); % switch mri order
opt.ana = opt.mri{1}.dat;
opt.clim = opt.mri{1}.clim; % use other scan's clim
set(opt.axes(4), 'Value', opt.clim(1)); % update minslider
set(opt.axes(5), 'Value', opt.clim(2)); % update maxslider
opt.slim = opt.mri{1}.slim; % use other scan's scatter clim
opt.axes(1:3) = opt.mri{1}.axes;
opt.init = true;

setappdata(h, 'opt', opt);
cb_redraw(h);
