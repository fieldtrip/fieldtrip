function [realign, snap] = ft_volumerealign(cfg, mri, target)

% FT_VOLUMEREALIGN spatially aligns an anatomical MRI with head coordinates
% based on external fiducials or anatomical landmarks. This function does
% not change the volume itself, but adjusts the homogeneous transformation
% matrix that describes the coordinate system. It also appends a
% coordsys-field to the mri, which specifies the coordinate system.
%
% This function only changes the coordinate system of an anatomical MRI, it
% does not change the MRI as such. For spatial normalisation (i.e. warping)
% of an MRI to a template brain you should use the FT_VOLUMENORMALISE
% function.
%
% Use as
%   [mri] = ft_volumerealign(cfg, mri)
%   [mri] = ft_volumerealign(cfg, mri, target);
% where the input mri should be a single anatomical or functional MRI
% volume that was for example read with FT_READ_MRI.
%
% The configuration can contain the following options
%   cfg.method         = different methods for aligning the volume
%                        'interactive', 'fiducial', 'landmark', 'headshape'
%                        'fsl', 'spm' (see below)
%   cfg.coordsys       = 'ctf' (default when specifying cfg.method =
%                         'interactive' or 'fiducial') or 'spm' (default
%                         when specifying cfg.method = 'landmark').
%                         Specifies the output coordinate system of the head. This
%                         string specifies the origin and the axes of the
%                         coordinate system. supported coordinate systems
%                         are: 'ctf', '4d', 'yokogawa', 'neuromag', 'itab'
%                         'spm', 'tal'.
%   cfg.clim           = [min max], scaling of the anatomy color (default
%                        is to adjust to the minimum and maximum)
%   cfg.parameter      = 'anatomy' the parameter which is used for the
%                         visualization
%
% When cfg.method = 'fiducial', the following cfg-option is required:
%   cfg.fiducial.nas  = [i j k], position of nasion
%   cfg.fiducial.lpa  = [i j k], position of LPA
%   cfg.fiducial.rpa  = [i j k], position of RPA
%   cfg.fiducial.zpoint = [i j k], a point on the positive z-axis. This is
%     an optional 'fiducial', and can be used to determine whether the
%     input voxel coordinate axes are left-handed (i.e. flipped in one of
%     the dimensions). If this additional point is specified, and the voxel
%     coordinate axes are left handed, the volume is flipped to yield right
%     handed voxel axes.
%
% When cfg.method = 'landmark', the following cfg-option is required:
%   cfg.landmark.ac      = [i j k], position of anterior commissure
%   cfg.landmark.pc      = [i j k], position of posterior commissure
%   cfg.landmark.xzpoint = [i j k], point on the midsagittal-plane with
%     positive Z-coordinate, i.e. interhemispheric point above ac and pc
% The coordinate system will be according to the RAS_Tal convention i.e.
% the origin corresponds with the anterior commissure the Y-axis is along
% the line from the posterior commissure to the anterior commissure the
% Z-axis is towards the vertex, in between the hemispheres the X-axis is
% orthogonal to the YZ-plane, positive to the right
%
% When cfg.method = 'interactive', a user interface allows for the
% specification of the fiducials or landmarks using the mouse, cursor keys
% and keyboard. Using the n/l/r keys the fiducials can be specified, the
% landmarks can be specified with a/p/z. When pressing q the interactive
% mode will stop and the transformation matrix is computed. This method
% also supports the cfg-option:
%  cfg.snapshot = 'no' ('yes'), making a snapshot of the image once a
%    fiducial or landmark location is selected. The optional second
%    output argument to the function will contain the handles to these
%    figures.
%  cfg.snapshotfile = 'ft_volumerealign_snapshot' or string, the root of
%    the filename for the snapshots, including the path. If no path is
%    given the files are saved to the pwd. The consecutive figures will be
%    numbered and saved as png-file.
%
% When cfg.method = 'headshape', the following cfg-option is required:
%  cfg.headshape = string pointing to a file describing a headshape, that
%    can be loaded with ft_read_headshape, or a fieldtrip-structure describing
%    a headshape
% 
% When cfg.method = 'fsl', a third input argument is required. The input volume is
% coregistered to this target volume, using fsl's flirt program. This
% assumes fsl to be installed. Options pertaining to the behavior of fsl
% should be defined in the subcfg cfg.fsl:
%   cfg.fsl.path    = string, specifying the path to fsl
%   cfg.fsl.costfun = string, specifying the cost-function used for
%                     coregistration 
%   cfg.fsl.interpmethod = string, specifying the interpolation method
%                     ('trilinear', 'nearestneighbour', 'sinc')
%   cfg.fsl.dof     = scalar, specifying the number of parameters for the
%                     affine transformation. 6 (rigid body), 7 (global
%                     rescale), 9 (traditional) or 12.
%   cfg.fsl.reslice = string, specifying whether the output image will be
%                     resliced conform the target image (default = 'yes') 
%
% When cfg.method = 'spm', a third input argument is required. The input volume is
% coregistered to this target volume, using spm. Options pertaining to the 
% behavior of spm can be defined in the subcfg cfg.spm:
%   cfg.spm.regtype = 'subj', 'rigid'
%   cfg.spm.smosrc  = scalar value 
%   cfg.spm.smoref  = scalar value
%
% With the 'interactive' and 'fiducial' methods it is possible to define an
% additional point (with the key 'z'), which should be a point on the
% positive side of the xy-plane, i.e. with a positive z-coordinate in world
% coordinates. This point will subsequently be used to check whether the
% input coordinate system is left or right-handed. For the 'interactive'
% and 'landmark' methods you can also specify an additional control point
% (with the key 'r'), that should be a point with a positive coordinate on
% the left-right axis.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ... cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a
% *.mat file on disk and/or the output data will be written to a *.mat
% file. These mat files should contain only a single variable,
% corresponding with the input/output structure.
%
% See also FT_READ_MRI, FT_ELECTRODEREALIGN, HEADCOORDINATES, SPM_AFFREG,
% SPM_NORMALISE

% Undocumented option:
%   cfg.weights = vector of weights that is used to weight the individual
%   headshape points in the icp algorithm. Used optionally in cfg.method     
%   = 'headshape'. If not specified, weights are put on points with
%   z-coordinate<0 (assuming those to be eye rims and nose ridges, i.e.
%   important points.
%
% Copyright (C) 2006-2011, Robert Oostenveld, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see
% http://www.ru.nl/neuroimaging/fieldtrip for the documentation and
% details.
%
%    FieldTrip is free software: you can redistribute it and/or modify it
%    under the terms of the GNU General Public License as published by the
%    Free Software Foundation, either version 3 of the License, or (at your
%    option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful, but
%    WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    General Public License for more details.
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
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar mri

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'realignfiducial', 'fiducial'});

% set the defaults
cfg.coordsys   = ft_getopt(cfg, 'coordsys',  '');
cfg.method     = ft_getopt(cfg, 'method',    ''); % deal with this below
cfg.fiducial   = ft_getopt(cfg, 'fiducial',  []);
cfg.landmark   = ft_getopt(cfg, 'landmark',  []);
cfg.parameter  = ft_getopt(cfg, 'parameter', 'anatomy');
cfg.clim       = ft_getopt(cfg, 'clim',      []);
cfg.snapshot   = ft_getopt(cfg, 'snapshot',  false);
cfg.snapshotfile = ft_getopt(cfg, 'snapshotfile', fullfile(pwd,'ft_volumerealign_snapshot'));

if strcmp(cfg.method, '')
  if isempty(cfg.landmark) && isempty(cfg.fiducial)
    cfg.method = 'interactive';
  elseif ~isempty(cfg.fiducial)
    cfg.method = 'fiducial';
  elseif ~isempty(cfg.landmark)
    cfg.method = 'landmark';
  end
end

if strcmp(cfg.coordsys, '')
  if strcmp(cfg.method, 'landmark')
    cfg.coordsys = 'spm';
  elseif strcmp(cfg.method, 'fiducial')
    cfg.coordsys = 'ctf';
  else
    cfg.coordsys = '';
  end
end

basedonmrk = strcmp(cfg.method, 'landmark');
basedonfid = strcmp(cfg.method, 'fiducial');

% these two have to be simultaneously true for a snapshot to be taken
dosnapshot   = istrue(cfg.snapshot);
if dosnapshot,
  % create an empty array of handles
  snap = [];
end
takesnapshot = false;

% select the parameter that should be displayed
cfg.parameter = parameterselection(cfg.parameter, mri);
if iscell(cfg.parameter) && ~isempty(cfg.parameter)
  cfg.parameter = cfg.parameter{1};
elseif iscell(cfg.parameter) && isempty(cfg.parameter)
  % cfg.parameter has been cleared by parameterselection due to a
  % dimensionality mismatch. Most probable cause here is the fact that a 4D
  % volume (e.g. DTI data) is in the input. This needs to be patched in a
  % more structural way at some point, but for the time being we'll use a
  % workaround here.
  
  % assume anatomy to be the parameter of interest
  siz = size(mri.anatomy);
  if all(siz(1:3)==mri.dim) && numel(siz)==4,
    % it's OK
    cfg.parameter= 'anatomy';
  else
    error('there''s an unexpected dimension mismatch');
  end
end

transform = [];
coordsys  = [];
switch cfg.method
  case 'fiducial'
    % do nothing
    
  case 'landmark'
    % do nothing
    
  case 'interactive'
    h  = figure;
    h1 = subplot('position',[0.02 0.55 0.44 0.44]);%subplot(2,2,1);
    h2 = subplot('position',[0.52 0.55 0.44 0.44]);%subplot(2,2,2);
    h3 = subplot('position',[0.02 0.05 0.44 0.44]);%subplot(2,2,3);
    handles = {h1 h2 h3};
    
    showcrosshair = true;
    showmarkers   = 1;
    dat = mri.(cfg.parameter);
    nas = [];
    lpa = [];
    rpa = [];
    antcomm = [];
    pstcomm = [];
    xzpoint = [];
    x = 1:mri.dim(1);
    y = 1:mri.dim(2);
    z = 1:mri.dim(3);
    xc = round(mri.dim(1)/2);
    yc = round(mri.dim(2)/2);
    zc = round(mri.dim(3)/2);
    
    updatepanel = [1 2 3];
    pnt         = zeros(0,3);
    markerpos   = zeros(0,3);
    markerlabel = {};
    markercolor = {};
    
    fprintf(strcat(...
      '1. To change the slice viewed in one plane, either:\n',...
      '   a. click (left mouse) in the image on a different plane. Eg, to view a more\n',...
      '      superior slice in the horizontal plane, click on a superior position in the\n',...
      '      coronal plane, or\n',...
      '   b. use the arrow keys to increase or decrease the slice number by one\n',...
      '2. To mark a fiducial position or anatomical landmark, do BOTH:\n',...
      '   (this can be done multiple times, until you are satisfied with the positions.\n',...
      '    for each type of point, the most recent selection is stored.)\n',...
      '   a. select the position by clicking on it in any slice with the left mouse\n',...
      '      button\n',...
      '   b. identify it by pressing either n/l/r for fiducials, or a/p/z for\n',...
      '      anatomical landmarks\n',...
      '   c. additional control point for the fiducials can be a point along the\n',...
      '      positive z-axis, press z\n',...
      '   d. additional control point for the landmarks can be a point along the\n',...
      '      positive x-axis (to the participant''s right), press r\n',...
      '3. To change the display:\n',...
      '   a. press c or C on keyboard to show/hide crosshair\n',...
      '   b. press m or M on keyboard to show/hide marked positions\n',...
      '   c. press + or - on (numeric) keyboard to change the color range''s upper limit\n',...
      '4. To finalize markers and quit interactive mode, press q on keyboard\n'));
    
    %'3. To unmark or remark a location\n',...
    %'   a. click with the middle mouse button to unmark last position\n',...
    %'   b. select new position with right mouse button and identify it using the\n',...
    %'      keyboard\n',...
    while true % break when 'q' is pressed
      %       fprintf('click with the left mouse button to reslice the display to a new position\n');
      %       fprintf('click with the right mouse button to mark a position\n');
      %       fprintf('click with the middle mouse button to unmark the last marked position\n');
      %       fprintf('press n/l/r on keyboard to record the current position as fiducial location\n');
      %       fprintf('press a/p/z on keyboard to record the current position as anatomical landmark\n');
      %       fprintf('press the arrow keys on the keyboard to increment or decrement the slice number by one\n');
      %       fprintf('press c or C on the keyboard to show or hide the crosshair\n');
      %       fprintf('press m or M on the keyboard to show or hide the marked positions\n');
      %       fprintf('press q on keyboard to quit interactive mode\n');
      
      
      xc = round(xc); xc = max(1,xc); xc = min(mri.dim(1),xc);
      yc = round(yc); yc = max(1,yc); yc = min(mri.dim(2),yc);
      zc = round(zc); zc = max(1,zc); zc = min(mri.dim(3),zc);
      markers = {markerpos markerlabel markercolor};
      [h1, h2, h3] = volplot(x, y, z, dat, [xc yc zc], cfg.clim, showcrosshair, updatepanel, handles, showmarkers, markers);
      drawnow;
      if dosnapshot && takesnapshot
        % create a new figure and draw right away, this will keep the old one on the screen
        snap(end+1) = copyobj(h(1), 0);
        set(snap(end), 'visible', 'off');
        print(snap(end), '-dpng', [cfg.snapshotfile,num2str(numel(snap))]);
        set(0, 'currentfigure', h(1));
      end
      takesnapshot = false;
      
      try, [d1, d2, key] = ginput(1); catch, key='q'; end
      switch key
        
        % contrast scaling
        case 43 % numpad +
          if isempty(cfg.clim)
            cfg.clim = [min(dat(:)) max(dat(:))];
          end
          % reduce color scale range by 10%
          cscalefactor = (cfg.clim(2)-cfg.clim(1))/10;
          cfg.clim(2) = cfg.clim(2)-cscalefactor;
        case 45 % numpad -
          if isempty(cfg.clim)
            cfg.clim = [min(dat(:)) max(dat(:))];
          end
          % increase color scale range by 10%
          cscalefactor = (cfg.clim(2)-cfg.clim(1))/10;
          cfg.clim(2) = cfg.clim(2)+cscalefactor;
          
        case 113 % 'q'
          delete(h(1));
          
          break;
        case 108 % 'l'
          lpa = [xc yc zc];
          takesnapshot = true;
        case 114 % 'r'
          rpa = [xc yc zc];
          takesnapshot = true;
        case 110 % 'n'
          nas = [xc yc zc];
          takesnapshot = true;
        case 97  % 'a'
          antcomm = [xc yc zc];
          takesnapshot = true;
        case 112 % 'p'
          pstcomm = [xc yc zc];
          takesnapshot = true;
        case 122 % 'z'
          xzpoint = [xc yc zc];
          takesnapshot = true;
        case 99  % 'c'
          showcrosshair = true;
        case 67  % 'C'
          showcrosshair = false;
        case 109 % 'm'
          showmarkers = 2;
        case 77 % 'M'
          showmarkers = 0;
        case 1 % left mouse click
          % update the view to a new position
          l1 = get(get(gca, 'xlabel'), 'string');
          l2 = get(get(gca, 'ylabel'), 'string');
          switch l1,
            case 'i'
              xc = d1;
            case 'j'
              yc = d1;
            case 'k'
              zc = d1;
            otherwise
              continue;
          end
          switch l2,
            case 'i'
              xc = d2;
            case 'j'
              yc = d2;
            case 'k'
              zc = d2;
            otherwise
              continue;
          end
          if l1=='i' && l2=='j'
            updatepanel = [1 2 3];
          elseif l1=='i' && l2=='k'
            updatepanel = [2 3 1];
          elseif l1=='j' && l2=='k'
            updatepanel = [3 1 2];
          end
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
          if l1=='i' && l2=='j'
            updatepanel = [1 2 3];
          elseif l1=='i' && l2=='k'
            updatepanel = [2 3 1];
          elseif l1=='j' && l2=='k'
            updatepanel = [3 1 2];
          end
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
        case 28 % arrow left
          % update the coordinates
          l1 = get(get(gca, 'xlabel'), 'string');
          l2 = get(get(gca, 'ylabel'), 'string');
          if l1=='i' && l2=='j'
            xc = xc-1; updatepanel = [1 2 3];
          elseif l1=='i' && l2=='k'
            xc = xc-1; updatepanel = [2 3 1];
          elseif l1=='j' && l2=='k'
            yc = yc-1; updatepanel = [3 1 2];
          end
        case 30 % arrow up
          % update the coordinates
          l1 = get(get(gca, 'xlabel'), 'string');
          l2 = get(get(gca, 'ylabel'), 'string');
          if l1=='i' && l2=='j'
            yc = yc+1; updatepanel = [1 2 3];
          elseif l1=='i' && l2=='k'
            zc = zc+1; updatepanel = [2 3 1];
          elseif l1=='j' && l2=='k'
            zc = zc+1; updatepanel = [3 1 2];
          end
        case 29 % arrow right
          % update the coordinates
          l1 = get(get(gca, 'xlabel'), 'string');
          l2 = get(get(gca, 'ylabel'), 'string');
          if l1=='i' && l2=='j'
            xc = xc+1; updatepanel = [1 2 3];
          elseif l1=='i' && l2=='k'
            xc = xc+1; updatepanel = [2 3 1];
          elseif l1=='j' && l2=='k'
            yc = yc+1; updatepanel = [3 1 2];
          end
        case 31 % arrow down
          % update the coordinates
          l1 = get(get(gca, 'xlabel'), 'string');
          l2 = get(get(gca, 'ylabel'), 'string');
          if l1=='i' && l2=='j'
            yc = yc-1; updatepanel = [1 2 3];
          elseif l1=='i' && l2=='k'
            zc = zc-1; updatepanel = [2 3 1];
          elseif l1=='j' && l2=='k'
            zc = zc-1; updatepanel = [3 1 2];
          end
        otherwise
          % do nothing
      end
      
      if all(round([xc yc zc])<=mri.dim) && all(round([xc yc zc])>0)
        fprintf('============================================================================\n');
        str = sprintf('voxel %d, indices [%d %d %d]', sub2ind(mri.dim(1:3), round(xc), round(yc), round(zc)), round([xc yc zc]));
        
        if isfield(mri, 'coordsys') && isfield(mri, 'unit')
          str = sprintf('%s, %s coordinates [%.1f %.1f %.1f] %s', str, mri.coordsys, warp_apply(mri.transform, [xc yc zc]), mri.unit);
        elseif ~isfield(mri, 'coordsys') && isfield(mri, 'unit')
          str = sprintf('%s, location [%.1f %.1f %.1f] %s', str, warp_apply(mri.transform, [xc yc zc]), mri.unit);
        elseif isfield(mri, 'coordsys') && ~isfield(mri, 'unit')
          str = sprintf('%s, %s coordinates [%.1f %.1f %.1f]', str, mri.coordsys, warp_apply(mri.transform, [xc yc zc]));
        elseif ~isfield(mri, 'coordsys') && ~isfield(mri, 'unis')
          str = sprintf('%s, location [%.1f %.1f %.1f]', str, warp_apply(mri.transform, [xc yc zc]));
        end
        fprintf('%s\n', str);
        % fprintf('cur_voxel = [%f %f %f], cur_head = [%f %f %f]\n', [xc yc zc], warp_apply(mri.transform, [xc yc zc]));
      end
      
      markerpos   = zeros(0,3);
      markerlabel = {};
      markercolor = {};
      if ~isempty(nas),
        fprintf('nas_voxel = [%f %f %f], nas_head = [%f %f %f]\n', nas, warp_apply(mri.transform, nas));
        markerpos   = [markerpos; nas];
        markerlabel = [markerlabel; {'nas'}];
        markercolor = [markercolor; {'b'}];
      end
      if ~isempty(lpa),
        fprintf('lpa_voxel = [%f %f %f], lpa_head = [%f %f %f]\n', lpa, warp_apply(mri.transform, lpa));
        markerpos   = [markerpos; lpa];
        markerlabel = [markerlabel; {'lpa'}];
        markercolor = [markercolor; {'g'}];
      end
      if ~isempty(rpa),
        fprintf('rpa_voxel = [%f %f %f], rpa_head = [%f %f %f]\n', rpa, warp_apply(mri.transform, rpa));
        markerpos   = [markerpos; rpa];
        markerlabel = [markerlabel; {'rpa'}];
        markercolor = [markercolor; {'r'}];
      end
      if ~isempty(antcomm),
        fprintf('antcomm_voxel = [%f %f %f], antcomm_head = [%f %f %f]\n', antcomm, warp_apply(mri.transform, antcomm));
        markerpos   = [markerpos; antcomm];
        markerlabel = [markerlabel; {'antcomm'}];
        markercolor = [markercolor; {'b'}];
      end
      if ~isempty(pstcomm),
        fprintf('pstcomm_voxel = [%f %f %f], pstcomm_head = [%f %f %f]\n', pstcomm, warp_apply(mri.transform, pstcomm));
        markerpos   = [markerpos; pstcomm];
        markerlabel = [markerlabel; {'pstcomm'}];
        markercolor = [markercolor; {'g'}];
      end
      if ~isempty(xzpoint),
        fprintf('xzpoint_voxel = [%f %f %f], xzpoint_head = [%f %f %f]\n', xzpoint, warp_apply(mri.transform, xzpoint));
        markerpos   = [markerpos; xzpoint];
        markerlabel = [markerlabel; {'xzpoint'}];
        markercolor = [markercolor; {'y'}];
      end
      if ~isempty(pnt)
        fprintf('%f extra points selected\n', size(pnt,1));
        markerpos   = [markerpos; pnt];
        markerlabel = [markerlabel; repmat({''}, size(pnt,1), 1)];
        markercolor = [markercolor; repmat({'m'}, size(pnt,1), 1)];
      end
    end % while true
    
    cfg.fiducial.nas    = nas;
    cfg.fiducial.lpa    = lpa;
    cfg.fiducial.rpa    = rpa;
    cfg.fiducial.zpoint = xzpoint;
    
    cfg.landmark.ac     = antcomm;
    cfg.landmark.pc     = pstcomm;
    cfg.landmark.xzpoint = xzpoint;
    cfg.landmark.rpoint  = rpa;
    
    if ~isempty(nas) && ~isempty(lpa) && ~isempty(rpa)
      basedonfid = 1;
      if isempty(cfg.coordsys)
        cfg.coordsys = 'ctf';
      end
    end
    
    if ~isempty(antcomm) && ~isempty(pstcomm) && ~isempty(xzpoint)
      basedonmrk = 1;
    end
    
  case 'headshape'
    if ischar(cfg.headshape)
      shape = ft_read_headshape(cfg.headshape);
    else
      shape = cfg.headshape;
    end
    shape = ft_convert_units(shape, 'mm');
    
    % extract skull surface from image
    tmpcfg        = [];
    tmpcfg.output = 'scalp';
    tmpcfg.smooth = 2;
    if isfield(cfg, 'template')
     tmpcfg.template = cfg.template;
    end
    seg           = ft_volumesegment(tmpcfg, mri);
    
    tmpcfg             = [];
    tmpcfg.method      = 'singleshell';
    tmpcfg.numvertices = 20000;
    scalp           = ft_prepare_headmodel(tmpcfg, seg);
    scalp           = ft_convert_units(scalp, 'mm');
    
    if ~isfield(cfg, 'weights')
      % weight the points with z-coordinate more than the rest. These are the
      % likely points that belong to the nose and eye rims
      w = ones(size(shape.pnt,1),1);
      %w(shape.pnt(:,3)<0) = 100; % this value seems to work
    else
      w = cfg.weights(:);
      if numel(w)~=size(shape.pnt,1),
        error('number of weights should be equal to the number of points in the headshape');
      end
    end
    
    % the icp function wants this as a function handle.
    weights = @(x)assignweights(x,w);
    
    % construct the coregistration matrix
    % [R, t, corr, D, data2] = icp2(scalp.bnd.pnt', shape.pnt', 20, [], weights); % icp2 is a helper function implementing the iterative closest point algorithm
    nrm         = normals(scalp.bnd.pnt, scalp.bnd.tri, 'vertex');
    [R, t, err] = icp(scalp.bnd.pnt', shape.pnt', 50, 'Minimize', 'plane', 'Normals', nrm', 'ReturnAll', true, 'Weight', weights, 'Extrapolation', true, 'WorstRejection', 0.05);
    [m,ix]      = min(err);
    R           = R(:,:,ix);
    t           = t(:,:,ix);
    transform   = inv([R t;0 0 0 1]);
    
    % warp the extracted scalp points to the new positions
    scalp.bnd.pnt          = warp_apply(transform, scalp.bnd.pnt);
    
    % create headshape structure for mri-based surface point cloud
    if isfield(mri, 'coordsys')
      scalp.coordsys = mri.coordsys;
    end
    
    % coordsys is the same as input mri
    coordsys = mri.coordsys;
    
    %     mrifid.pnt   = warp_apply(transform*mri.transform, [fiducials.nas;fiducials.lpa;fiducials.rpa]);
    %     mrifid.label = {'NZinteractive';'Linteractive';'Rinteractive'};
    %     shapemri.fid = mrifid;
       
    % update the cfg
    cfg.headshape    = shape;
    cfg.headshapemri = scalp;
    
    % touch it to survive trackconfig
    cfg.headshapemri;
    
  case 'fsl'
    if ~isfield(cfg, 'fsl'), cfg.fsl = []; end
    cfg.fsl.path         = ft_getopt(cfg.fsl, 'path',    '');
    cfg.fsl.costfun      = ft_getopt(cfg.fsl, 'costfun', 'corratio');
    cfg.fsl.interpmethod = ft_getopt(cfg.fsl, 'interpmethod', 'trilinear');
    cfg.fsl.dof          = ft_getopt(cfg.fsl, 'dof',     6);
    cfg.fsl.reslice      = ft_getopt(cfg.fsl, 'reslice', 'yes');
    cfg.fsl.searchrange  = ft_getopt(cfg.fsl, 'searchrange', [-180 180]);
    
    % write the input and target to a temporary file
    % and create some additional temporary file names to contain the output
    tmpname1 = tempname;
    tmpname2 = tempname;
    tmpname3 = tempname;
    tmpname4 = tempname;
    
    tmpcfg = [];
    tmpcfg.parameter = 'anatomy';
    tmpcfg.filename  = tmpname1;
    tmpcfg.filetype  = 'nifti';
    fprintf('writing the input volume to a temporary file: %s\n', [tmpname1,'.nii']);
    ft_volumewrite(tmpcfg, mri);
    tmpcfg.filename  = tmpname2;
    fprintf('writing the  target volume to a temporary file: %s\n', [tmpname2,'.nii']);
    ft_volumewrite(tmpcfg, target);
    
    % create the command to call flirt
    fprintf('using flirt for the coregistration\n');
    r1  = num2str(cfg.fsl.searchrange(1));
    r2  = num2str(cfg.fsl.searchrange(2));
    str = sprintf('%s/flirt -in %s -ref %s -out %s -omat %s -bins 256 -cost %s -searchrx %s %s -searchry %s %s -searchrz %s %s -dof %s -interp %s',...
          cfg.fsl.path, tmpname1, tmpname2, tmpname3, tmpname4, cfg.fsl.costfun, r1, r2, r1, r2, r1, r2, num2str(cfg.fsl.dof), cfg.fsl.interpmethod);
    if isempty(cfg.fsl.path), str = str(2:end); end % remove the first filesep, assume path to flirt to be known
        
    % system call
    system(str);
    
    % process the output
    if ~istrue(cfg.fsl.reslice)
      % get the transformation that corresponds to the coregistration and
      % reconstruct the mapping from the target's world coordinate system
      % to the input's voxel coordinate system
      
      fid = fopen(tmpname4);
      tmp = textscan(fid, '%f');
      fclose(fid);
      
      % this transforms from input voxels to target voxels
      vox2vox = reshape(tmp{1},4,4)';
      
      if det(target.transform(1:3,1:3))>0
        % flirt apparently flips along the x-dim if the det < 0
        % if images are not radiological, the x-axis is flipped, see:
        %  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0810&L=FSL&P=185638
        %  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0903&L=FSL&P=R93775

        % flip back
        flipmat = eye(4); flipmat(1,1) = -1; flipmat(1,4) = target.dim(1); 
        vox2vox = flipmat*vox2vox;
      end
      if det(mri.transform(1:3,1:3))>0
        % flirt apparently flips along the x-dim if the det < 0
        % flip back
        flipmat = eye(4); flipmat(1,1) = -1; flipmat(1,4) = mri.dim(1); 
        vox2vox = vox2vox*flipmat;
      end
      
      % very not sure about this (e.g. is vox2vox really doing what I think
      % it is doing? should I care about 0 and 1 based conventions?)
      % changing handedness?
      mri.transform = target.transform*vox2vox;
      
      transform = eye(4);
      if isfield(target, 'coordsys')
        coordsys = target.coordsys;
      else
        coordsys = 'unknown';
      end
      
    else
      % get the updated anatomy
      mrinew        = ft_read_mri([tmpname3, '.nii.gz']);
      mri.anatomy   = mrinew.anatomy;
      mri.transform = mrinew.transform;
      mri.dim       = mrinew.dim;
      
      transform = eye(4);
      if isfield(target, 'coordsys')
        coordsys = target.coordsys;
      else
        coordsys = 'unknown';
      end
    end
    delete([tmpname1,'.nii']);
    delete([tmpname2,'.nii']);
    delete([tmpname3,'.nii.gz']);
    delete(tmpname4);
    
  case 'spm'
    % ensure spm8 on the path
    ft_hastoolbox('SPM8', 1);

    if ~isfield(cfg, 'spm'), cfg.spm = []; end
    cfg.spm.regtype = ft_getopt(cfg.spm, 'regtype', 'subj');
    cfg.spm.smosrc  = ft_getopt(cfg.spm, 'smosrc',  2);
    cfg.spm.smoref  = ft_getopt(cfg.spm, 'smoref',  2);
    
    if ~isfield(mri,    'coordsys'), 
      mri = ft_convert_coordsys(mri); 
    else
      fprintf('Input volume has coordinate system ''%s''\n', mri.coordsys);
    end
    if ~isfield(target, 'coordsys'),
      target = ft_convert_coordsys(target);
    else
      fprintf('Target volume has coordinate system ''%s''\n', target.coordsys);
    end
    if strcmp(mri.coordsys, target.coordsys)
      % this should hopefully work
    else
      % only works when it is possible to approximately align the input to
      % the target coordsys
      if strcmp(target.coordsys, 'spm')
        mri = ft_convert_coordsys(mri, 'spm');
      else
        error('The coordinate systems of the input and target volumes are different, coregistration is not possible');
      end
    end
    
    % flip and permute the 3D volume itself, so that the voxel and
    % headcoordinates approximately correspond
    [tmp,    pvec_mri,    flip_mri, T] = align_ijk2xyz(mri);
    [target]                           = align_ijk2xyz(target);
        
    tname1 = [tempname, '.img'];
    tname2 = [tempname, '.img'];
    V1 = ft_write_mri(tname1, tmp.anatomy,    'transform', tmp.transform,    'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
    V2 = ft_write_mri(tname2, target.anatomy, 'transform', target.transform, 'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  
    flags         = cfg.spm;
    flags.nits    = 0; %set number of non-linear iterations to zero
    params        = spm_normalise(V2,V1,[],[],[],flags);
    %mri.transform = (target.transform/params.Affine)/T;
    transform = (target.transform/params.Affine)/T/mri.transform;
    % transform     = eye(4);
    if isfield(target, 'coordsys')
      coordsys = target.coordsys;
    else
      coordsys = 'unknown';
    end
    
    % delete the temporary files
    delete(tname1);
    delete(tname2);
  
  otherwise
    error('unsupported method');
end

if basedonfid && basedonmrk
  basedonmrk = 0;
  warning('both fiducials and anatomical landmarks have been defined interactively: using the fiducials for realignment');
end

if basedonfid
  % the fiducial locations are now specified in voxels, convert them to head
  % coordinates according to the existing transform matrix
  nas_head = warp_apply(mri.transform, cfg.fiducial.nas);
  lpa_head = warp_apply(mri.transform, cfg.fiducial.lpa);
  rpa_head = warp_apply(mri.transform, cfg.fiducial.rpa);
  if isfield(cfg.fiducial, 'zpoint') && ~isempty(cfg.fiducial.zpoint)
    zpnt_head = warp_apply(mri.transform, cfg.fiducial.zpoint);
    [transform, coordsys] = headcoordinates(nas_head, lpa_head, rpa_head, zpnt_head, cfg.coordsys);
  else
    % compute the homogeneous transformation matrix describing the new coordinate system
    [transform, coordsys] = headcoordinates(nas_head, lpa_head, rpa_head, cfg.coordsys);
  end
elseif basedonmrk
  % the fiducial locations are now specified in voxels, convert them to head
  % coordinates according to the existing transform matrix
  ac     = warp_apply(mri.transform, cfg.landmark.ac);
  pc     = warp_apply(mri.transform, cfg.landmark.pc);
  xzpoint= warp_apply(mri.transform, cfg.landmark.xzpoint);
  if isfield(cfg.landmark, 'rpoint') && ~isempty(cfg.landmark.rpoint)
    rpnt_head = warp_apply(mri.transform, cfg.landmark.rpoint);
    [transform, coordsys] = headcoordinates(ac, pc, xzpoint, rpnt_head, 'spm');
  else
    % compute the homogenous transformation matrix describing the new coordinate system
    [transform, coordsys] = headcoordinates(ac, pc, xzpoint, 'spm');
  end
  
else
  % something else has created a transform and coordsys
  
end % if basedonXXX

% copy the input anatomical or functional volume
realign = mri;

if ~isempty(transform)
  % combine the additional transformation with the original one
  realign.transformorig = mri.transform;
  realign.transform     = transform * mri.transform;
  realign.coordsys      = coordsys;
else
  warning('no coordinate system realignment has been done');
end

if exist('pnt', 'var') && ~isempty(pnt)
  realign.marker = pnt;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous mri
ft_postamble history realign
ft_postamble savevar realign


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to show three orthogonal slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h1, h2, h3] = volplot(x, y, z, dat, c, cscale, showcrosshair, updatepanel, handles, showmarkers, markers)

xi = c(1);
yi = c(2);
zi = c(3);

% manual color scaling of anatomy data is usefull in case of some pixel noise
if nargin<6 || isempty(cscale)
  cmin = min(dat(:));
  cmax = max(dat(:));
else
  cmin = cscale(1);
  cmax = cscale(2);
end

if nargin<8
  updatepanel = [1 2 3];
end

if nargin<9
  h1 = [];
  h2 = [];
  h3 = [];
else
  h1 = handles{1};
  h2 = handles{2};
  h3 = handles{3};
end

if showmarkers==1
  markerpos = round(markers{1});
  markercolor = markers{3};
  sel1 = find(markerpos(:,2)==repmat(c(2),size(markerpos,1),1));
  sel2 = find(markerpos(:,1)==repmat(c(1),size(markerpos,1),1));
  sel3 = find(markerpos(:,3)==repmat(c(3),size(markerpos,1),1));
elseif showmarkers==2
  markerpos = round(markers{1});
  markercolor = markers{3};
  sel1 = 1:size(markerpos,1);
  sel2 = 1:size(markerpos,1);
  sel3 = 1:size(markerpos,1);
end

for k = 1:numel(updatepanel)
  
  update = updatepanel(k);
  if update==1
    subplot(h1);
    imagesc(x, z, squeeze(dat(:,yi,:))'); set(gca, 'ydir', 'normal')
    xlabel('i'); ylabel('k');
    caxis([cmin cmax]);
    if showcrosshair
      crosshair([x(xi) z(zi)], 'color', 'yellow');
    end
    if showmarkers && numel(sel1)>0
      hold on;
      for kk = 1:numel(sel1)
        plot(markerpos(sel1(kk),1), markerpos(sel1(kk),3), 'marker', '.', 'color', markercolor{sel1(kk)});
      end
      hold off;
    end
    axis equal; axis tight;
  end
  
  if update==2
    subplot(h2);
    imagesc(y, z, squeeze(dat(xi,:,:))'); set(gca, 'ydir', 'normal')
    xlabel('j'); ylabel('k');
    caxis([cmin cmax]);
    if showcrosshair
      crosshair([y(yi) z(zi)], 'color', 'yellow');
    end
    if showmarkers && numel(sel2)>0
      hold on;
      for kk = 1:numel(sel2)
        plot(markerpos(sel2(kk),2), markerpos(sel2(kk),3), 'marker', '.', 'color', markercolor{sel2(kk)});
      end
      hold off;
    end
    axis equal; axis tight;
  end
  
  if update==3
    subplot(h3);
    imagesc(x, y, squeeze(dat(:,:,zi))'); set(gca, 'ydir', 'normal')
    xlabel('i'); ylabel('j');
    caxis([cmin cmax]);
    if showcrosshair
      crosshair([x(xi) y(yi)], 'color', 'yellow');
    end
    if showmarkers && numel(sel3)>0
      hold on;
      for kk = 1:numel(sel3)
        plot(markerpos(sel3(kk),1), markerpos(sel3(kk),2), 'marker', '.', 'color', markercolor{sel3(kk)});
      end
      hold off;
    end
    axis equal; axis tight;
  end
  
end

colormap gray

h = gca;

function [R, t, corr, error, data2] = icp2(data1, data2, res, tri, weights)

% [R, t, corr, error, data2] = icp2(data1, data2, res, tri)
%
% This is an implementation of the Iterative Closest Point (ICP) algorithm.
% The function takes two data sets and registers data2 with data1. It is
% assumed that data1 and data2 are in approximation registration. The code
% iterates till no more correspondences can be found.
%
% This is a modified version (12 April, 2005). It is more accurate and has
% less chances of getting stuck in a local minimum as opposed to my earlier
% version icp.m
%
% Arguments: data1 - 3 x n matrix of the x, y and z coordinates of data set 1
%            data2 - 3 x m matrix of the x, y and z coordinates of data set 2
%            res   - the tolerance distance for establishing closest point
%                     correspondences. Normally set equal to the resolution
%                     of data1
%            tri   - optional argument. obtained by tri = delaunayn(data1');
%
% Returns: R - 3 x 3 accumulative rotation matrix used to register data2
%          t - 3 x 1 accumulative translation vector used to register data2
%          corr - p x 3 matrix of the index no.s of the corresponding points of
%                 data1 and data2 and their corresponding Euclidean distance
%          error - the mean error between the corresponding points of data1
%                  and data2 (normalized with res)
%          data2 - 3 x m matrix of the registered data2
%
%
% Copyright : This code is written by Ajmal Saeed Mian {ajmal@csse.uwa.edu.au}
%              Computer Science, The University of Western Australia. The code
%              may be used, modified and distributed for research purposes with
%              acknowledgement of the author and inclusion of this copyright information.

maxIter = 500;
c1 = 0;
c2 = 1;
R = eye(3);
t = zeros(3,1);
if nargin < 4 || isempty(tri)
    tri = delaunayn(data1');
end
n = 0;
while c2 ~= c1
  c1 = c2;
  [corr, D] = dsearchn(data1', tri, data2');
  corr(:,2:3)     = [(1 : length(corr))' D];
  corr(D>2*res,:) = [];
  
  corr = -sortrows(-corr,3);
  corr = sortrows(corr,1);
  [B, Bi, Bj] = unique(corr(:,1));
  corr = corr(Bi,:);
  
  [R1, t1] = reg(data1, data2, corr, weights);
  data2 = R1*data2;
  data2 = [data2(1,:)+t1(1); data2(2,:)+t1(2); data2(3,:)+t1(3)];
  R = R1*R;
  t = R1*t + t1;
  c2 = length(corr);
  n = n + 1;
  if n > maxIter
    break;
  end
end

e1 = 1000001;
e2 = 1000000;
n = 0;
noChangeCount = 0;
while noChangeCount < 10
  e1 = e2;
  [corr, D] = dsearchn(data1', tri, data2');
  corr(:,2:3) = [(1:length(corr))' D];
  corr(D>2*res,:) = [];
  
  corr = -sortrows(-corr,3);
  corr = sortrows(corr,1);
  [B, Bi, Bj] = unique(corr(:,1));
  corr = corr(Bi,:);
  
  [R1 t1] = reg(data1, data2, corr, weights);
  data2 = R1*data2;
  data2 = [data2(1,:)+t1(1); data2(2,:)+t1(2); data2(3,:)+t1(3)];
  R = R1*R;
  t = R1*t + t1;
  e2 = sum(corr(:,3))/(length(corr)*res);
  
  n = n + 1;
  if n > maxIter
    break;
  end
  if abs(e2-e1)<res/10000
    noChangeCount = noChangeCount + 1;
  end
end
error = min(e1,e2);

%-----------------------------------------------------------------
function [R1, t1] = reg(data1, data2, corr, weights)

n = length(corr);
if nargin<4
  weights = ones(n,1);
end
M = data1(:,corr(:,1));
mm = mean(M,2);
S = data2(:,corr(:,2));%*sparse(diag(weights(corr(:,2))));
ms = mean(S,2);
Sshifted = [S(1,:)-ms(1); S(2,:)-ms(2); S(3,:)-ms(3)];
Mshifted = [M(1,:)-mm(1); M(2,:)-mm(2); M(3,:)-mm(3)];
K = Sshifted*sparse(diag(weights(corr(:,2))))*Mshifted';
K = K/n;
[U A V] = svd(K);
R1 = V*U';
if det(R1)<0
  B = eye(3);
  B(3,3) = det(V*U');
  R1 = V*B*U';
end
t1 = mm - R1*ms;


%-------------------------------------------------------
function y = assignweights(x, w)

% x is an indexing vector with the same number of arguments as w
y = w(:)';
