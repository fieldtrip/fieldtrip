function [elec] = ft_electrodeplacement(cfg, varargin)

% FT_ELECTRODEPLACEMENT allows placing electrodes on volume or headshape.
% The different methods are described in detail below.

% VOLUME - Navigate an orthographic display of a volume (e.g. CT or
% MR scan), and assign an electrode label to the current crosshair location
% by clicking on a label in the eletrode list. You can undo the selection by 
% clicking on the same label again. The electrode labels shown in the list 
% can be prespecified using cfg.channel when calling ft_electrodeplacement.
% The zoom slider allows zooming in at the location of the crosshair.
% The intensity sliders allow thresholding the image's low and high values.
% The magnet feature transports the crosshair to the nearest peak intensity 
% voxel, within a 3 voxel radius of the selected location. 
% The labels feature displays the labels of the selected electrodes within 
% the orthoplot.
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

% The configuration can contain the following options
%   cfg.method         = string representing the method for aligning or placing the electrodes
%                        'mri'             place electrodes in a brain volume
%                        'headshape'       place electrodes on the head surface
%   cfg.channel        = Nx1 cell-array with selection of channels (default = '1','2', ...)
%
% See also FT_ELECTRODEREALIGN, FT_VOLUMEREALIGN

% Copyright (C) 2015, Arjen Stolk & Robert Oostenveld
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


% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar mri
ft_preamble provenance mri
ft_preamble trackconfig

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% set the defaults
cfg.parameter  = ft_getopt(cfg, 'parameter', 'anatomy');
cfg.method     = ft_getopt(cfg, 'method');               % volume, headshape
cfg.channel    = ft_getopt(cfg, 'channel',          []); % default will be determined further down {'1', '2', ...}
cfg.clim       = ft_getopt(cfg, 'clim',          [0 1]); % initial volume intensity limit
cfg.radius     = ft_getopt(cfg, 'radius',            2); % magnet feature: detect peaks within n radius physical voxels   
cfg.elec       = ft_getopt(cfg, 'elec',             []); % already defined electrodes, for plotting purposes 

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
  case 'headshape'
    headshape = fixpos(varargin{1});
    headshape = ft_determine_coordsys(headshape);
end

switch cfg.method
  case 'headshape'
    % give the user instructions
    disp('Use the mouse to click on the desired electrode positions');
    disp('Afterwards you may have to update the electrode labels');
    disp('Press "q" when you are done');
    % open a figure
    figure;
    % plot the faces of the 2D or 3D triangulation
    skin = [255 213 119]/255;
    ft_plot_mesh(headshape,'facecolor', skin,'EdgeColor','none','facealpha',0.7);
    lighting gouraud
    material shiny
    camlight
    % rotate3d on
    xyz = ft_select_point3d(headshape, 'nearest', false, 'multiple', true, 'marker', '*');
    numelec = size(xyz, 1);
    
    % construct the output electrode structure
    elec = keepfields(headshape, {'unit', 'coordsys'});
    elec.elecpos = xyz;
    for i=1:numelec
      try
        elec.label{i} = cfg.channel{i};
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
    
    % add callbacks
    set(h, 'windowbuttondownfcn', @cb_buttonpress);
    set(h, 'windowbuttonupfcn',   @cb_buttonrelease);
    set(h, 'windowkeypressfcn',   @cb_keyboard);
    set(h, 'CloseRequestFcn',     @cb_cleanup);
    
    % axes settings
    xdim = mri.dim(1) + mri.dim(2);
    ydim = mri.dim(2) + mri.dim(3);
    
    xsize(1) = 0.82*mri.dim(1)/xdim;
    xsize(2) = 0.82*mri.dim(2)/xdim;
    ysize(1) = 0.82*mri.dim(3)/ydim;
    ysize(2) = 0.82*mri.dim(2)/ydim;
    
    xc = round(mri.dim(1)/2); % start view
    yc = round(mri.dim(2)/2);
    zc = round(mri.dim(3)/2);
    
    % axis handles will hold the anatomical data if present, along with labels etc.
    h1 = axes('position',[0.07 0.07+ysize(2)+0.05 xsize(1) ysize(1)]);
    h2 = axes('position',[0.07+xsize(1)+0.05 0.07+ysize(2)+0.05 xsize(2) ysize(1)]);
    h3 = axes('position',[0.07 0.07 xsize(1) ysize(2)]);
    set(h1,'Tag','ik','Visible','off','XAxisLocation','top');
    set(h2,'Tag','jk','Visible','off','XAxisLocation','top');
    set(h3,'Tag','ij','Visible','off');
    
    % inspect transform matrix, if the voxels are isotropic then the screen
    % pixels also should be square
    hasIsotropicVoxels = norm(mri.transform(1:3,1)) == norm(mri.transform(1:3,2))...
      && norm(mri.transform(1:3,2)) == norm(mri.transform(1:3,3));
    if hasIsotropicVoxels
      set(h1,'DataAspectRatio',[1 1 1]);
      set(h2,'DataAspectRatio',[1 1 1]);
      set(h3,'DataAspectRatio',[1 1 1]);
    end
    
    dat = double(mri.(cfg.parameter));
    dmin = min(dat(:));
    dmax = max(dat(:));
    dat  = (dat-dmin)./(dmax-dmin); % range between 0 and 1
    
    % intensity range sliders
    h45text = uicontrol('Style', 'text',...
      'String','Intensity',...
      'Units', 'normalized', ...
      'Position',[2*xsize(1)+0.03 ysize(2)+0.03 xsize(1)/4 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on');
    
    h4 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 1, ...
      'Value', cfg.clim(1), ...
      'Units', 'normalized', ...
      'Position', [2*xsize(1)+0.02 0.10+ysize(2)/3 0.05 ysize(2)/2], ...
      'Callback', @cb_minslider);
    
    h5 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 1, ...
      'Value', cfg.clim(2), ...
      'Units', 'normalized', ...
      'Position', [2*xsize(1)+0.07 0.10+ysize(2)/3 0.05 ysize(2)/2], ...
      'Callback', @cb_maxslider);
    
    % intensity range slider (dual-knob slider): the java component gives issues when wanting to
    % access the opt structure
    % [jRangeSlider] = com.jidesoft.swing.RangeSlider(0,1,cfg.clim(1),cfg.clim(2));  % min,max,low,high
    % [jRangeSlider, h4] = javacomponent(jRangeSlider, [], h);
    % set(h4, 'Units', 'normalized', 'Position', [0.05+xsize(1) 0.07 0.07 ysize(2)], 'Parent', h);
    % set(jRangeSlider, 'Orientation', 1, 'PaintTicks', true, 'PaintLabels', true, ...
    %     'Background', java.awt.Color.white, 'StateChangedCallback', @cb_intensityslider);
    
    % slice slider
    h6text = uicontrol('Style', 'text',...
      'String','Slice',...
      'Units', 'normalized', ...
      'Position',[0.07 0.02 xsize(1)/3-0.03 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on');
    
    h6 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 1, 'Max', mri.dim(3), ...
      'Value', zc, ...
      'Units', 'normalized', ...
      'Position', [xsize(1)/3+0.05 0.01 2*(xsize(1)/3) 0.05], ...
      'SliderStep', [1/mri.dim(3) 1/mri.dim(3)], ...
      'Callback', @cb_sliceslider);
    
    % electrode listbox
    if isempty(cfg.channel) % ability to specify grid types, e.g. grid1_elec1, strip1_elec1?
      for c = 1:150
        cfg.channel{c} = sprintf('%d', c);
      end
    end
    for c = 1:numel(cfg.channel)
      chanstrings{c} = ['<HTML><FONT color="silver">' cfg.channel{c} '</FONT></HTML>']; % hmtl'ize
    end
    h7 = uicontrol('Style', 'listbox', ...
      'Parent', h, ...
      'Value', [], 'Min', 0, 'Max', numel(chanstrings), ...
      'Units', 'normalized', ...
      'Position', [0.07+xsize(1)+0.05 0.07 xsize(1)/2 ysize(2)], ...
      'Callback', @cb_eleclistbox, ...
      'String', chanstrings);
    
    % switches / radio buttons
    h8 = uicontrol('Style', 'radiobutton',...
      'Parent', h, ...
      'Value', 1, ...
      'String','Magnet',...
      'Units', 'normalized', ...
      'Position',[2*xsize(1) 0.15 xsize(1)/3 0.05],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_magnetbutton);
    
    h9 = uicontrol('Style', 'radiobutton',...
      'Parent', h, ...
      'Value', 0, ...
      'String','Labels',...
      'Units', 'normalized', ...
      'Position',[2*xsize(1) 0.07 xsize(1)/3 0.05],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on', ...
      'Callback', @cb_labelsbutton);
    
    % intensity range sliders
    h10text = uicontrol('Style', 'text',...
      'String','Zoom',...
      'Units', 'normalized', ...
      'Position',[1.8*xsize(1)+0.01 ysize(2)+0.03 xsize(1)/4 0.04],...
      'BackgroundColor', [1 1 1], ...
      'HandleVisibility','on');
    
    h10 = uicontrol('Style', 'slider', ...
      'Parent', h, ...
      'Min', 0, 'Max', 0.9, ...
      'Value', 0, ...
      'Units', 'normalized', ...
      'Position', [1.8*xsize(1)+0.02 0.10+ysize(2)/3 0.05 ysize(2)/2], ...
      'SliderStep', [.1 .1], ...
      'Callback', @cb_zoomslider);
    
    if ~isempty(cfg.elec)
      for e = 1:numel(cfg.elec.label)
        markers{e,1} = cfg.elec.voxorig(e,:);
        markers{e,2} = cfg.elec.elecpos(e,:);
        markers{e,3} = cfg.elec.label{e};
      end
    else
      markervox   = zeros(0,3);
      markerpos   = zeros(0,3);
      markerlabel = {};
      markers = repmat({markervox markerpos markerlabel},numel(cfg.channel),1);
    end
    
    % instructions to the user
    fprintf(strcat(...
      ' \n'));
    
    % create structure to be passed to gui
    opt               = [];
    opt.dim           = mri.dim;
    opt.ijk           = [xc yc zc];
    opt.xsize         = xsize;
    opt.ysize         = ysize;
    opt.handlesaxes   = [h1 h2 h3 h4 h5 h6 h7 h8 h9 h10];
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
    opt.magnet        = get(h8, 'Value');
    opt.radius        = cfg.radius;
    opt.showmarkers   = true;
    opt.markers       = markers;
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
    elec.voxorig = [];
    for i=1:length(opt.markers)
      if ~isempty(opt.markers{i,1})
        elec.label = [elec.label; opt.markers{i,3}];
        elec.elecpos = [elec.elecpos; opt.markers{i,2}];
        elec.voxorig = [elec.voxorig; opt.markers{i,1}]; % keep the original voxel coordinates
      end
    end
    elec.chanpos  = elec.elecpos; % identicial to elecpos
    elec.tra = eye(size(elec.elecpos,1));
    if isfield(mri, 'coordsys')
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
xyz = mri.transform * opt.ijk;
opt.ijk = opt.ijk(1:3)';

% construct a string with user feedback
str1 = sprintf('voxel %d, index [%d %d %d]', sub2ind(mri.dim(1:3), xi, yi, zi), opt.ijk);

if opt.init
  ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', opt.ijk, 'style', 'subplot', 'parents', [h1 h2 h3], 'update', opt.update, 'doscale', false, 'clim', opt.clim);
  
  opt.anahandles = findobj(opt.handlesfigure, 'type', 'surface')';
  parenttag = get(opt.anahandles,'parent');
  parenttag{1} = get(parenttag{1}, 'tag');
  parenttag{2} = get(parenttag{2}, 'tag');
  parenttag{3} = get(parenttag{3}, 'tag');
  [i1,i2,i3] = intersect(parenttag, {'ik';'jk';'ij'});
  opt.anahandles = opt.anahandles(i3(i2)); % seems like swapping the order
  opt.anahandles = opt.anahandles(:)';
  set(opt.anahandles, 'tag', 'ana');
  
  opt.h1axis = axis(h1); % for zooming purposes
  opt.h2axis = axis(h2);
  opt.h3axis = axis(h3);
else
  ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', opt.ijk, 'style', 'subplot', 'surfhandle', opt.anahandles, 'update', opt.update, 'doscale', false, 'clim', opt.clim);
  
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
xloadj = round((xi-opt.h1axis(1))-(xi-opt.h1axis(1))*opt.zoom)+.5;
xhiadj = round((opt.h1axis(2)-xi)-(opt.h1axis(2)-xi)*opt.zoom)+.5;
yloadj = round((yi-opt.h1axis(3))-(yi-opt.h1axis(3))*opt.zoom)+.5;
yhiadj = round((opt.h1axis(4)-yi)-(opt.h1axis(4)-yi)*opt.zoom)+.5;
zloadj = round((zi-opt.h1axis(5))-(zi-opt.h1axis(5))*opt.zoom)+.5;
zhiadj = round((opt.h1axis(6)-zi)-(opt.h1axis(6)-zi)*opt.zoom)+.5;
axis(h1, [xi-xloadj xi+xhiadj yi-yloadj yi+yhiadj zi-zloadj zi+zhiadj]);
axis(h2, [xi-xloadj xi+xhiadj yi-yloadj yi+yhiadj zi-zloadj zi+zhiadj]);
axis(h3, [xi-xloadj xi+xhiadj yi-yloadj yi+yhiadj]);

if opt.init
  % draw the crosshairs for the first time
  hch1 = crosshair([xi yi zi], 'parent', h1, 'color', 'yellow'); % [xi 1 zi]
  hch2 = crosshair([xi yi zi], 'parent', h2, 'color', 'yellow'); % [opt.dim(1) yi zi]
  hch3 = crosshair([xi yi zi], 'parent', h3, 'color', 'yellow'); % [xi yi opt.dim(3)]
  opt.handlescross  = [hch1(:)';hch2(:)';hch3(:)'];
  opt.handlesmarker = [];
else
  % update the existing crosshairs, don't change the handles
  crosshair([xi yi zi], 'handle', opt.handlescross(1, :));
  crosshair([xi yi zi], 'handle', opt.handlescross(2, :));
  crosshair([xi yi zi], 'handle', opt.handlescross(3, :));
end

if opt.showcrosshair
  set(opt.handlescross,'Visible','on');
else
  set(opt.handlescross,'Visible','off');
end

delete(opt.handlesmarker(opt.handlesmarker(:)>0));
opt.handlesmarker = [];

for i=1:size(opt.markers,1)
  if ~isempty(opt.markers{i,1})
    pos = opt.markers{i,1}; % voxel coordinates
    
    posi = pos(1);
    posj = pos(2);
    posk = pos(3);
    
    subplot(h1);
    hold on
    opt.handlesmarker(i,1) = plot3(posi, posj, posk, 'marker', '+', 'color', 'r'); % posi, 1, posk
    if opt.showlabels
      opt.handlesmarker(i,4) = text(posi, posj, posk, opt.markers{i,3});
    end
    hold off
    
    subplot(h2);
    hold on
    opt.handlesmarker(i,2) = plot3(posi, posj, posk, 'marker', '+', 'color', 'r'); % opt.dim(1), posj, posk
    if opt.showlabels
      opt.handlesmarker(i,5) = text(posi, posj, posk, opt.markers{i,3});
    end
    hold off
    
    subplot(h3);
    hold on
    opt.handlesmarker(i,3) = plot3(posi, posj, posk, 'marker', '+', 'color', 'r'); % posi, posj, opt.dim(3)
    if opt.showlabels
      opt.handlesmarker(i,6) = text(posi, posj, posk, opt.markers{i,3});
    end
    hold off
  end
end % for each marker
  
% adjust slice slider accordingly
set(opt.handlesaxes(6), 'Value', opt.ijk(3));

opt.init = false;
setappdata(h, 'opt', opt);
set(h, 'currentaxes', curr_ax);

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
    
  case {'i' 'j' 'k' 'm' 28 29 30 31 'leftarrow' 'rightarrow' 'uparrow' 'downarrow'} % TODO FIXME use leftarrow rightarrow uparrow downarrow
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

h   = getparent(h);
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

h   = getparent(h);
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
    radius = opt.radius;
    cubic = opt.ana(center(1)-radius:center(1)+radius, center(2)-radius:center(2)+radius, center(3)-radius:center(3)+radius);
    [val, idx] = max(cubic(:)); % find peak intensity voxel within the radius
    [ix, iy, iz] = ind2sub(size(cubic), idx);
    opt.ijk = center+[ix, iy, iz]-radius-1;
  end
end
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_cleanup(h, eventdata)

opt = getappdata(h, 'opt');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_minslider(h4, eventdata)

newlim = get(h4, 'value');
h = getparent(h4);
opt = getappdata(h, 'opt');
opt.clim(1) = newlim;
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_maxslider(h5, eventdata)

newlim = get(h5, 'value');
h = getparent(h5);
opt = getappdata(h, 'opt');
opt.clim(2) = newlim;
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_intensityslider(h4, eventdata)

loval = get(h4, 'value');
hival = get(h4, 'highvalue');
h = getparent(h4); % this fails: The name 'parent' is not an accessible property for an instance of class 'com.jidesoft.swing.RangeSlider'.
opt = getappdata(h, 'opt');
opt.clim = [loval hival];
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_sliceslider(h6, eventdata)

newslice = round(get(h6, 'value'));
h = getparent(h6);
opt = getappdata(h, 'opt');
opt.ijk(3) = newslice;
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_eleclistbox(h7, eventdata)

elecidx = get(h7, 'Value'); % chosen elec
if ~isempty(elecidx)
  eleclis = cellstr(get(h7, 'String')); % all labels
  eleclab = eleclis{elecidx}; % this elec's label
  
  h = getparent(h7);
  opt = getappdata(h, 'opt');
  
  % toggle electrode status and assign markers
  if strfind(eleclab, 'silver') % not yet, check
    eleclab = regexprep(eleclab, '"silver"','"black"'); % replace font color
    opt.markers(elecidx,:) = {opt.vox opt.pos opt.label(elecidx)};  % assign marker position and label
  elseif strfind(eleclab, 'black') % already chosen before, uncheck
    eleclab = regexprep(eleclab, '"black"','"silver"'); % replace font color
    opt.markers(elecidx,:) = {zeros(0,3) zeros(0,3) {}};  % assign marker position and label
  end
  
  % update plot
  eleclis{elecidx} = eleclab;
  set(h7, 'String', eleclis);
  setappdata(h, 'opt', opt);
  cb_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_magnetbutton(h8, eventdata)

h = getparent(h8);
opt = getappdata(h, 'opt');
opt.magnet = get(h8, 'value');
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_labelsbutton(h9, eventdata)

h = getparent(h9);
opt = getappdata(h, 'opt');
opt.showlabels = get(h9, 'value');
setappdata(h, 'opt', opt);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_zoomslider(h10, eventdata)

h = getparent(h10);
opt = getappdata(h, 'opt');
opt.zoom = round(get(h10, 'value')*10)/10;
setappdata(h, 'opt', opt);
cb_redraw(h);