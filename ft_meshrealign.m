function [mesh_realigned] = ft_meshrealign(cfg, mesh)

% FT_MESHREALIGN rotates, translates and optionally scales electrode positions. The
% different methods are described in detail below.
%
% INTERACTIVE - You can display the mesh surface together with axis coordinate
% system, and manually (using the graphical user interface) adjust the rotation,
% translation and scaling parameters.
%
% FIDUCIAL - The coordinate system is updated according to the definition of the
% coordinates of anatomical landmarks or fiducials that are specified in the
% configuration. If the fiducials are not specified in the configurartion, you will
% have to click them in an interactive display of the mesh surface.
%
% Use as
%   mesh = ft_meshrealign(cfg, mesh)
% where the mesh input argument comes from FT_READ_HEADSHAPE or FT_PREPARE_MESH and
% cfg is a configuration structure that should contain
%
%  cfg.method    = string, can be 'interactive' or fiducial' (default = 'interactive')
%
% The configuration can furthermore contain
%   cfg.coordsys        = string, can be 'ctf', 'neuromag', '4d', 'bti', 'itab'
%   cfg.fiducial.nas    = [x y z], position of nasion
%   cfg.fiducial.lpa    = [x y z], position of LPA
%   cfg.fiducial.rpa    = [x y z], position of RPA
%
% The fiducials should be expressed in the coordinates and units of the input mesh.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_READ_HEADSHAPE, FT_PREPARE_MESH, FT_ELECTRODEREALIGN, FT_VOLUMEREALIGN

% Copyrights (C) 2017, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    mesh
ft_preamble provenance mesh
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
mesh = ft_checkdata(mesh, 'datatype', 'mesh', 'feedback', 'yes', 'hasunit', 'yes', 'hascoordsys', 'no');

% get the options
cfg.method    = ft_getopt(cfg, 'method', 'interactive');
cfg.coordsys  = ft_getopt(cfg, 'coordsys');
cfg.fiducial = ft_getopt(cfg, 'fiducial');
cfg.fiducial.nas = ft_getopt(cfg.fiducial, 'nas');
cfg.fiducial.lpa = ft_getopt(cfg.fiducial, 'lpa');
cfg.fiducial.rpa = ft_getopt(cfg.fiducial, 'rpa');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start with a copy
mesh_realigned = keepfields(mesh, {'pos', 'tri', 'tet', 'hex', 'unit', 'line', 'edge', 'color', 'curv', 'sulc'});

switch cfg.method
  case 'interactive'
    
    tmpcfg = [];
    tmpcfg.template.axes = 'yes';
    tmpcfg.template.headshape.pos       = zeros(3,3); % three vertices
    tmpcfg.template.headshape.tri       = [1 2 3];    % one triangle
    tmpcfg.template.headshape.unit      = 'mm';
    tmpcfg.template.headshape.coordsys  = cfg.coordsys;
    tmpcfg.template.headshapestyle = {'vertexcolor', 'none', 'edgecolor', 'none', 'facecolor', 'none'};
    tmpcfg.individual.headshape = mesh_realigned;
    tmpcfg = ft_interactiverealign(tmpcfg);
    % keep the homogenous transformation
    transform = tmpcfg.m;
    
  case 'fiducial'
    % it must be present and it must be a string
    ft_checkopt(cfg, 'coordsys', {'char'});
    
    hasnas = ~isempty(cfg.fiducial.nas);
    haslpa = ~isempty(cfg.fiducial.lpa);
    hasrpa = ~isempty(cfg.fiducial.rpa);
    
    if ~hasnas || ~haslpa || ~hasrpa
      % do something interactive to get them
      % this is shared with ft_volumerealign
      
      cfg.fiducial.nas    = [nan nan nan];
      cfg.fiducial.lpa    = [nan nan nan];
      cfg.fiducial.rpa    = [nan nan nan];
      cfg.fiducial.zpoint = [nan nan nan];
      
      scalp = mesh;
      mri = [];
      
      switch cfg.coordsys
        case {'ctf' '4d' 'bti' 'yokogawa' 'asa' 'itab' 'neuromag'}
          fidlabel  = {'nas', 'lpa', 'rpa', 'zpoint'};
          fidletter = {'n', 'l', 'r', 'z'};
          fidexplanation1 = '      press n for nas, l for lpa, r for rpa\n';
          fidexplanation2 = '      press z for an extra control point that should have a positive z-value\n';
        case 'acpc'
          fidlabel  = {'ac', 'pc', 'xzpoint', 'right'};
          fidletter = {'a', 'p', 'z', 'r'};
          fidexplanation1 = '      press a for ac, p for pc, z for xzpoint\n';
          fidexplanation2 = '      press r for an extra control point that should be on the right side\n';
        case 'paxinos'
          fidlabel  = {'bregma', 'lambda', 'yzpoint'};
          fidletter = {'b', 'l', 'z'};
          fidexplanation1 = '      press b for bregma, l for lambda, z for yzpoint\n';
          fidexplanation2 = '';
        otherwise
          ft_error('unknown coordinate system "%s"', cfg.coordsys);
      end % switch coordsys
      
      fprintf('\n');
      fprintf(strcat(...
        '1. To change the orientation of the head surface, use the\n',...
        '"Rotate 3D" option in the figure toolbar\n',...
        '2. To mark a fiducial position or anatomical landmark, do BOTH:\n',...
        '   a. select the position by clicking on it with the left mouse button\n',...
        '   b. specify it by pressing the letter corresponding to the fiducial/landmark:\n', fidexplanation1, fidexplanation2, ...
        '   You can mark the fiducials multiple times, until you are satisfied with the positions.\n',...
        '3. To finalize markers and quit interactive mode, press q on keyboard\n'));
      
      % start building the figure
      h = figure;
      set(h, 'color', [1 1 1]);
      set(h, 'visible', 'on');
      % add callbacks
      set(h, 'windowkeypressfcn',   @cb_keyboard_surface);
      set(h, 'CloseRequestFcn',     @cb_quit);
      
      % create figure handles
      h1 = axes;
      
      % create structure to be passed to gui
      opt                 = [];
      opt.viewresult      = false; % flag to use for certain keyboard/redraw calls
      opt.handlesfigure   = h;
      opt.handlesaxes     = h1;
      opt.handlesfigure   = h;
      opt.handlesmarker   = [];
      opt.camlighthandle  = [];
      opt.init            = true;
      opt.quit            = false;
      opt.scalp           = scalp;
      opt.showmarkers     = false;
      opt.mri             = mri;
      opt.fiducial        = cfg.fiducial;
      opt.fidlabel        = fidlabel;
      opt.fidletter       = fidletter;
      opt.fidexplanation1 = fidexplanation1;
      if isfield(scalp, 'unit') && ~strcmp(scalp.unit, 'unknown')
        opt.unit = scalp.unit;  % this is shown in the feedback on screen
      else
        opt.unit = '';        % this is not shown
      end
      
      setappdata(h, 'opt', opt);
      cb_redraw_surface(h);
      
      while(opt.quit==0)
        uiwait(h);
        opt = getappdata(h, 'opt');
      end
      delete(h);
      
      % store the interactively determined fiducials in the configuration
      % the actual coordinate transformation will be done further down
      cfg.fiducial = opt.fiducial;
      
    end
    
    % compute the homogenous transformation
    transform = ft_headcoordinates(cfg.fiducial.nas, cfg.fiducial.lpa, cfg.fiducial.rpa, cfg.coordsys);
    
    
  otherwise
    ft_error('unsupported method "%s"', cfg.method);
end

% update the positions
mesh_realigned = ft_transform_geometry(transform, mesh_realigned);

% assign the coordinate system
if ~isempty(cfg.coordsys)
  mesh_realigned.coordsys = cfg.coordsys;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   mesh
ft_postamble provenance mesh_realigned
ft_postamble history    mesh_realigned
ft_postamble savevar    mesh_realigned


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw_surface(h, eventdata)
h   = getparent(h);
opt = getappdata(h, 'opt');

markercolor = {'r', 'g', 'b', 'y'};

if opt.init
  ft_plot_mesh(opt.scalp, 'edgecolor', 'none', 'facecolor', 'skin')
  hold on
end

% recreate the camera lighting
delete(opt.camlighthandle);
opt.camlighthandle = camlight;

% remove the previous fiducials
delete(opt.handlesmarker(opt.handlesmarker(:)>0));
opt.handlesmarker = [];

% redraw the fiducials
for i=1:length(opt.fidlabel)
  lab = opt.fidlabel{i};
  pos = opt.fiducial.(lab);
  if all(~isnan(pos))
    opt.handlesmarker(i,1) = plot3(pos(1), pos(2), pos(3), 'marker', 'o', 'color', markercolor{i});
    opt.handlesmarker(i,2) = text(pos(1), pos(2), pos(3), lab);
  end
end

opt.init = false;
setappdata(h, 'opt', opt);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard_surface(h, eventdata)
h   = getparent(h);
opt = getappdata(h, 'opt');

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end

% get the most recent surface position that was clicked with the mouse
pos = select3d(opt.handlesaxes);

sel = find(strcmp(opt.fidletter, key));
if ~isempty(sel)
  % update the corresponding fiducial
  opt.fiducial.(opt.fidlabel{sel}) = pos(:)';
end

fprintf('==================================================================================\n');
for i=1:length(opt.fidlabel)
  lab = opt.fidlabel{i};
  vox = opt.fiducial.(lab);
  pos = vox;
  ind = nan;
  switch opt.unit
    case 'mm'
      fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.1f %.1f %.1f] %s\n', lab, ind, round(vox), pos, opt.unit);
    case 'cm'
      fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.2f %.2f %.2f] %s\n', lab, ind, round(vox), pos, opt.unit);
    case 'm'
      fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.4f %.4f %.4f] %s\n', lab, ind, round(vox), pos, opt.unit);
    otherwise
      fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%f %f %f] %s\n', lab, ind, round(vox), pos, opt.unit);
  end
end

setappdata(h, 'opt', opt);

if isequal(key, 'q')
  cb_quit(h);
else
  cb_redraw_surface(h);
end

uiresume(h);


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

% extract to-be-plotted/clicked location and check whether inside figure
xi = opt.ijk(1);
yi = opt.ijk(2);
zi = opt.ijk(3);
if any([xi yi zi] > mri.dim) || any([xi yi zi] <= 0)
  return;
end

% transform here to coordinate system space instead of voxel space if viewing results
% the code were this transform will impact fiducial/etc coordinates is unaffected, as it is switched off
% (note: fiducial/etc coordinates are transformed into coordinate space in the code dealing with realignment)
if opt.viewresult
  tmp = ft_warp_apply(mri.transform, [xi yi zi]);
  xi = tmp(1);
  yi = tmp(2);
  zi = tmp(3);
end

if opt.init
  % create the initial figure
  if ~opt.viewresult
    % if realigning, plotting is done in voxel space
    ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', [xi yi zi], 'style', 'subplot', 'parents', [h1 h2 h3], 'update', opt.update, 'doscale', false, 'clim', opt.clim);
  else
    % if viewing result, plotting is done in head coordinate system space
    if ~opt.twovol
      % one vol case
      ft_plot_ortho(opt.ana, 'transform', mri.transform, 'location', [xi yi zi], 'style', 'subplot', 'parents', [h1 h2 h3], 'update', opt.update, 'doscale', false, 'clim', opt.realignclim);
    else
      % two vol case
      % base volume, with color red
      hbase = []; % need the handle for the individual surfs
      [hbase(1), hbase(2), hbase(3)] = ft_plot_ortho(opt.ana, 'transform', mri.transform, 'unit', mri.unit, 'location', [xi yi zi], 'style', 'subplot', 'parents', [h1 h2 h3], 'update', opt.update, 'doscale', false, 'clim', opt.targetclim, 'datmask',opt.targetmask, 'opacitylim', [0 1]);
      for ih = 1:3
        col = get(hbase(ih), 'CData');
        col(:,:,2:3) = 0;
        set(hbase(ih), 'CData',col);
      end
      % aligned volume, with color blue
      hreal = []; % need the handle for the individual surfs
      [hreal(1), hreal(2), hreal(3)] = ft_plot_ortho(opt.realignana, 'transform', opt.realignvol.transform, 'unit', opt.realignvol.unit, 'location', [xi yi zi], 'style', 'subplot', 'parents', [h1 h2 h3], 'update', opt.update, 'doscale', false, 'clim', opt.realignclim, 'datmask',opt.realignmask, 'opacitylim', [0 1]);
      for ih = 1:3
        col = get(hreal(ih), 'CData');
        col(:,:,1:2) = 0;
        set(hreal(ih), 'CData',col);
      end
    end
  end % if ~opt.viewresult
  
  % fetch surf objects, set ana tag, and put in surfhandles
  if ~opt.viewresult || (opt.viewresult && ~opt.twovol)
    opt.anahandles = findobj(opt.handlesfigure, 'type', 'surface')';
    parenttag = get(opt.anahandles, 'parent');
    parenttag{1} = get(parenttag{1}, 'tag');
    parenttag{2} = get(parenttag{2}, 'tag');
    parenttag{3} = get(parenttag{3}, 'tag');
    [i1,i2,i3] = intersect(parenttag, {'ik';'jk';'ij'});
    opt.anahandles = opt.anahandles(i3(i2)); % seems like swapping the order
    opt.anahandles = opt.anahandles(:)';
    set(opt.anahandles, 'tag', 'ana');
  else
    % this should do the same as the above
    set(hbase, 'tag', 'ana');
    set(hreal, 'tag', 'ana');
    opt.anahandles = {hbase, hreal};
  end
else
  % update the existing figure
  if ~opt.viewresult
    % if realigning, plotting is done in voxel space
    ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', [xi yi zi], 'style', 'subplot', 'surfhandle', opt.anahandles, 'update', opt.update, 'doscale', false, 'clim', opt.clim);
  else
    % if viewing result, plotting is done in head coordinate system space
    if ~opt.twovol
      % one vol case
      ft_plot_ortho(opt.ana, 'transform', mri.transform, 'unit', mri.unit, 'location', [xi yi zi], 'style', 'subplot', 'surfhandle', opt.anahandles, 'update', opt.update, 'doscale', false, 'clim', opt.realignclim);
    else
      % two vol case
      % base volume, with color red
      hbase = []; % need the handle for the individual surfs
      [hbase(1), hbase(2), hbase(3)] = ft_plot_ortho(opt.ana, 'transform', mri.transform, 'unit', mri.unit, 'location', [xi yi zi], 'style', 'subplot', 'surfhandle', opt.anahandles{1}, 'update', opt.update, 'doscale', false, 'clim', opt.targetclim, 'datmask', opt.targetmask, 'opacitylim', [0 1]);
      for ih = 1:3
        col = get(hbase(ih), 'CData');
        col(:,:,2:3) = 0;
        set(hbase(ih), 'CData', col);
      end
      % aligned volume, with color blue
      hreal = []; % need the handle for the individual surfs
      [hreal(1), hreal(2), hreal(3)] = ft_plot_ortho(opt.realignana, 'transform', opt.realignvol.transform, 'unit', opt.realignvol.unit, 'location', [xi yi zi], 'style', 'subplot', 'surfhandle', opt.anahandles{2}, 'update', opt.update, 'doscale', false, 'clim', opt.realignclim, 'datmask', opt.realignmask, 'opacitylim', [0 1]);
      for ih = 1:3
        col = get(hreal(ih), 'CData');
        col(:,:,1:2) = 0;
        set(hreal(ih), 'CData', col);
      end
    end
  end % if ~opt.viewresult
  
  % display current location
  if ~opt.viewresult
    % if realigning, plotting is done in voxel space
    if all(round([xi yi zi])<=mri.dim) && all(round([xi yi zi])>0)
      fprintf('==================================================================================\n');
      
      lab = 'crosshair';
      vox = [xi yi zi];
      ind = sub2ind(mri.dim(1:3), round(vox(1)), round(vox(2)), round(vox(3)));
      pos = ft_warp_apply(mri.transform, vox);
      switch opt.unit
        case 'mm'
          fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.1f %.1f %.1f] %s\n', lab, ind, vox, pos, opt.unit);
        case 'cm'
          fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.2f %.2f %.2f] %s\n', lab, ind, vox, pos, opt.unit);
        case 'm'
          fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.4f %.4f %.4f] %s\n', lab, ind, vox, pos, opt.unit);
        otherwise
          fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%f %f %f] %s\n', lab, ind, vox, pos, opt.unit);
      end
    end
    
    for i=1:length(opt.fidlabel)
      lab = opt.fidlabel{i};
      vox = opt.fiducial.(lab);
      ind = sub2ind(mri.dim(1:3), round(vox(1)), round(vox(2)), round(vox(3)));
      pos = ft_warp_apply(mri.transform, vox);
      switch opt.unit
        case 'mm'
          fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.1f %.1f %.1f] %s\n', lab, ind, vox, pos, opt.unit);
        case 'cm'
          fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.2f %.2f %.2f] %s\n', lab, ind, vox, pos, opt.unit);
        case 'm'
          fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%.4f %.4f %.4f] %s\n', lab, ind, vox, pos, opt.unit);
        otherwise
          fprintf('%10s: voxel %9d, index = [%3d %3d %3d], head = [%f %f %f] %s\n', lab, ind, vox, pos, opt.unit);
      end
    end
    
  else
    % if viewing result, plotting is done in head coordinate system space
    lab = 'crosshair';
    pos = [xi yi zi];
    switch opt.unit
      case 'mm'
        fprintf('%10s: head = [%.1f %.1f %.1f] %s\n', lab, pos, opt.unit);
      case 'cm'
        fprintf('%10s: head = [%.2f %.2f %.2f] %s\n', lab, pos, opt.unit);
      case 'm'
        fprintf('%10s: head = [%.4f %.4f %.4f] %s\n', lab, pos, opt.unit);
      otherwise
        fprintf('%10s: head = [%f %f %f] %s\n', lab, pos, opt.unit);
    end
  end % if ~opt.viewresult
  
end % if opt.init

set(opt.handlesaxes(1), 'Visible', 'on');
set(opt.handlesaxes(2), 'Visible', 'on');
set(opt.handlesaxes(3), 'Visible', 'on');
if opt.viewresult
  set(opt.handlesaxes(1), 'color', [.94 .94 .94]);
  set(opt.handlesaxes(2), 'color', [.94 .94 .94]);
  set(opt.handlesaxes(3), 'color', [.94 .94 .94]);
end


% make the last current axes current again
sel = findobj('type', 'axes', 'tag',tag);
if ~isempty(sel)
  set(opt.handlesfigure, 'currentaxes', sel(1));
end

% set crosshair coordinates dependent on voxel/system coordinate space
% crosshair needs to be plotted 'towards' the viewing person, i.e. with a little offset
% i.e. this is the coordinate of the 'flat' axes with a little bit extra in the direction of the axis
% this offset cannot be higher than the to be plotted data, or it will not be visible (i.e. be outside of the visible axis)
if ~opt.viewresult
  crossoffs = opt.dim;
  crossoffs(2) = 1; % workaround to use the below
else
  % because the orientation of the three slices are determined by eye(3) (no orientation is specified above),
  % the direction of view is always:
  % h1 -to+
  % h2 +to-
  % h3 -to+
  % use this to create the offset for viewing the crosshair
  mincoordstep = abs(ft_warp_apply(mri.transform, [1 1 1]) - ft_warp_apply(mri.transform, [2 2 2]));
  crossoffs = [xi yi zi] + [1 -1 1].*mincoordstep;
end

if opt.init
  % draw the crosshairs for the first time
  hch1 = ft_plot_crosshair([xi crossoffs(2) zi], 'parent', h1, 'color', 'yellow');
  hch2 = ft_plot_crosshair([crossoffs(1) yi zi], 'parent', h2, 'color', 'yellow');
  hch3 = ft_plot_crosshair([xi yi crossoffs(3)], 'parent', h3, 'color', 'yellow');
  opt.handlescross  = [hch1(:)';hch2(:)';hch3(:)'];
  opt.handlesmarker = [];
else
  % update the existing crosshairs, don't change the handles
  ft_plot_crosshair([xi crossoffs(2) zi], 'handle', opt.handlescross(1, :));
  ft_plot_crosshair([crossoffs(1) yi zi], 'handle', opt.handlescross(2, :));
  ft_plot_crosshair([xi yi crossoffs(3)], 'handle', opt.handlescross(3, :));
end
% For some unknown god-awful reason, the line command 'disables' all transparency.
% The below command resets it. It was the only axes property that I (=roemei) could
% find that changed after adding the crosshair, and putting it back to 'childorder'
% instead of 'depth' fixes the problem. Lucky, the line command only 'disables' in
% the new graphics system introduced in 2014b (any version below is fine, and does
% not contain the sortmethod property --> crash)
if ~verLessThan('matlab', '8.4') % 8.4 = 2014b
  set(h1, 'sortMethod', 'childorder')
  set(h2, 'sortMethod', 'childorder')
  set(h3, 'sortMethod', 'childorder')
end

if opt.showcrosshair
  set(opt.handlescross, 'Visible', 'on');
else
  set(opt.handlescross, 'Visible', 'off');
end

markercolor = {'r', 'g', 'b', 'y'};

delete(opt.handlesmarker(opt.handlesmarker(:)>0));
opt.handlesmarker = [];

if ~opt.viewresult
  for i=1:length(opt.fidlabel)
    pos = opt.fiducial.(opt.fidlabel{i});
    %   if any(isnan(pos))
    %     continue
    %   end
    
    posi = pos(1);
    posj = pos(2);
    posk = pos(3);
    
    subplot(h1);
    hold on
    opt.handlesmarker(i,1) = plot3(posi, 1, posk, 'marker', 'o', 'color', markercolor{i});
    hold off
    
    subplot(h2);
    hold on
    opt.handlesmarker(i,2) = plot3(opt.dim(1), posj, posk, 'marker', 'o', 'color', markercolor{i});
    hold off
    
    subplot(h3);
    hold on
    opt.handlesmarker(i,3) = plot3(posi, posj, opt.dim(3), 'marker', 'o', 'color', markercolor{i});
    hold off
  end % for each fiducial
end

if opt.showmarkers
  set(opt.handlesmarker, 'Visible', 'on');
else
  set(opt.handlesmarker, 'Visible', 'off');
end

opt.init = false;
setappdata(h, 'opt', opt);
set(h, 'currentaxes', curr_ax);

uiresume


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
switch key
  case {'' 'shift+shift' 'alt-alt' 'control+control' 'command-0'}
    % do nothing
    
  case '1'
    subplot(opt.handlesaxes(1));
    
  case '2'
    subplot(opt.handlesaxes(2));
    
  case '3'
    subplot(opt.handlesaxes(3));
    
  case opt.fidletter
    if ~opt.viewresult
      sel = strcmp(key, opt.fidletter);
      fprintf('==================================================================================\n');
      fprintf('selected %s\n', opt.fidlabel{sel});
      opt.fiducial.(opt.fidlabel{sel}) = opt.ijk;
      setappdata(h, 'opt', opt);
      cb_redraw(h);
    end
    
  case 'q'
    setappdata(h, 'opt', opt);
    cb_quit(h);
    
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
    end
    
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
  case {43 'add' 'shift+equal'}  % + or numpad +
    % contrast scaling
    % disable if viewresult
    if ~opt.viewresult
      if isempty(opt.clim)
        opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
      end
      % reduce color scale range by 10%
      cscalefactor = (opt.clim(2)-opt.clim(1))/10;
      %opt.clim(1) = opt.clim(1)+cscalefactor;
      opt.clim(2) = opt.clim(2)-cscalefactor;
      setappdata(h, 'opt', opt);
      cb_redraw(h);
    end
    
  case {45 'subtract' 'hyphen' 'shift+hyphen'} % - or numpad -
    % contrast scaling
    % disable if viewresult
    if ~opt.viewresult
      if isempty(opt.clim)
        opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
      end
      % increase color scale range by 10%
      cscalefactor = (opt.clim(2)-opt.clim(1))/10;
      %opt.clim(1) = opt.clim(1)-cscalefactor;
      opt.clim(2) = opt.clim(2)+cscalefactor;
      setappdata(h, 'opt', opt);
      cb_redraw(h);
    end
    
  case 99  % 'c'
    opt.showcrosshair = ~opt.showcrosshair;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
  case 102 % 'f'
    if ~opt.viewresult
      opt.showmarkers = ~opt.showmarkers;
      setappdata(h, 'opt', opt);
      cb_redraw(h);
    end
    
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
    
  otherwise
    % do nothing
    
end % switch key
if ~opt.viewresult
  uiresume(h)
end
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

uiresume

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
uiresume


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_getposition(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h,       'currentaxes');
pos     = mean(get(curr_ax, 'currentpoint'));

tag = get(curr_ax, 'tag');

% transform pos from coordinate system space to voxel space if viewing results
if opt.viewresult
  pos = ft_warp_apply(inv(opt.mri.transform),pos); % not sure under which circumstances the transformation matrix is not invertible...
end

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

setappdata(h, 'opt', opt);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quit(h, eventdata)

opt = getappdata(h, 'opt');
if ~opt.viewresult
  opt.quit = true;
  setappdata(h, 'opt', opt);
  uiresume
else
  % not part of interactive process requiring output handling, quite immediately
  delete(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

