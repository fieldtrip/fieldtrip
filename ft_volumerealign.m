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
% or
%   [mri] = ft_volumerealign(cfg, mri, target)
% where the input mri should be a single anatomical or functional MRI volume that was
% for example read with FT_READ_MRI and the optional target MRI is the target
% anatomical MRI for FSL.
%
% The configuration can contain the following options
%   cfg.method         = different methods for aligning the volume
%                        'interactive', 'fiducial', 'landmark', 'headshape'
%                        'fsl', 'spm' (see below)
%   cfg.coordsys       = 'ctf' (default when specifying cfg.method =
%                         'interactive' or 'fiducial') or 'spm' (default
%                         when specifying cfg.method = 'landmark').
%                         Specifies the output coordinate system of the head.
%                         This string specifies the origin and the axes of the
%                         coordinate system. supported coordinate systems
%                         are 'ctf', '4d', 'yokogawa', 'neuromag', 'itab'
%                         'spm' and 'tal'.
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
% When cfg.method = 'headshape', the function extracts the scalp surface from
% the input MRI, and aligns this surface with the headshape. Options pertaining
% to this method can be defined in the subcfg cfg.headshape. The following
% cfg-option is required:
%  cfg.headshape.headshape = string pointing to a file describing a headshape, that
%    can be loaded with FT_READ_HEADSHAPE, or a FieldTrip-structure describing
%    a headshape
%
% The following options are optional:
%  cfg.headshape.scalpsmooth    = scalar (default = 2): smoothing parameter
%    for the scalp extraction
%  cfg.headsahpe.scalpthreshold = scalar (default = 0.1): threshold parameter
%    for the scalp extraction
%  cfg.headshape.interactive    = 'yes' ('no'): use interactive realignment
%    to align headshape with scalp surface
%  cfg.headshape.icp            = 'yes' ('no'): use automatic realignment
%    based on the icp-algorithm. If both 'interactive' and 'icp' are
%    executed, the icp step follows the interactive realignment step.
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

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'realignfiducial', 'fiducial'});
cfg = ft_checkconfig(cfg, 'renamed', {'landmark', 'fiducial'});
cfg = ft_checkconfig(cfg, 'renamedval', {'landmark', 'fiducial'});

% set the defaults
cfg.coordsys   = ft_getopt(cfg, 'coordsys',  []);
cfg.method     = ft_getopt(cfg, 'method',    []); % deal with this below
cfg.fiducial   = ft_getopt(cfg, 'fiducial',  []);
cfg.parameter  = ft_getopt(cfg, 'parameter', 'anatomy');
cfg.clim       = ft_getopt(cfg, 'clim',      []);
cfg.snapshot   = ft_getopt(cfg, 'snapshot',  false);
cfg.snapshotfile = ft_getopt(cfg, 'snapshotfile', fullfile(pwd,'ft_volumerealign_snapshot'));

if isempty(cfg.method)
  if isempty(cfg.fiducial)
    % fiducials have not yet been specified
    cfg.method = 'interactive';
  else
    % fiducials have already been specified
    cfg.method = 'fiducial';
  end
end

if isempty(cfg.coordsys)
  if     isstruct(cfg.fiducial) && all(ismember(fieldnames(cfg.fiducial), {'lpa', 'rpa', 'nas', 'zpoint'}))
    cfg.coordsys = 'ctf';
  elseif isstruct(cfg.fiducial) && all(ismember(fieldnames(cfg.fiducial), {'ac', 'pc', 'xzpoint', 'right'}))
    cfg.coordsys = 'spm';
  elseif strcmp(cfg.method, 'interactive')
    cfg.coordsys = 'ctf';
  end
  warning('defaulting to %s coordinate system', cfg.coordsys);
end

% these two have to be simultaneously true for a snapshot to be taken
dosnapshot = istrue(cfg.snapshot);
if dosnapshot,
  % create an empty array of handles
  snap = [];
end

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

% start with an empty transform and coordsys
transform = [];
coordsys  = [];

if strcmp(cfg.method, 'fiducial') || strcmp(cfg.method, 'interactive')
  switch cfg.coordsys
    case {'ctf' '4d' 'bti' 'yokogawa' 'asa' 'itab' 'neuromag'}
      fidlabel  = {'nas', 'lpa', 'rpa', 'zpoint'};
      fidletter = {'n', 'l', 'r', 'z'};
      fidexplanation1 = '      press n for nas, l for lpa, r for rpa\n';
      fidexplanation2 = '      press z for an extra control point that should have a positive z-value\n';
    case {'spm' 'tal'}
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
      error('unknown coordinate system "%s"', cfg.coordsys);
  end
  
  for i=1:length(fidlabel)
    if ~isfield(cfg.fiducial, fidlabel{i}) || isempty(cfg.fiducial.(fidlabel{i}))
      cfg.fiducial.(fidlabel{i}) = [nan nan nan];
    end
  end
end % interactive or fiducial

switch cfg.method
  case 'fiducial'
    % the actual coordinate transformation will be done further down
    
  case 'landmark'
    % the actual coordinate transformation will be done further down
    
  case 'interactive'
    %% start building the figure
    h = figure;
    set(h, 'color', [1 1 1]);
    set(h, 'visible', 'on');
    
    % add callbacks
    set(h, 'windowbuttondownfcn', @cb_buttonpress);
    set(h, 'windowbuttonupfcn',   @cb_buttonrelease);
    set(h, 'windowkeypressfcn',   @cb_keyboard);
    set(h, 'CloseRequestFcn',     @cb_cleanup);
    
    % enforce the size of the subplots to be isotropic
    xdim = mri.dim(1) + mri.dim(2);
    ydim = mri.dim(2) + mri.dim(3);
    
    % inspect transform matrix, if the voxels are isotropic then the screen
    % pixels also should be square
    hasIsotropicVoxels = norm(mri.transform(1:3,1)) == norm(mri.transform(1:3,2))...
      && norm(mri.transform(1:3,2)) == norm(mri.transform(1:3,3));
    
    xsize(1) = 0.82*mri.dim(1)/xdim;
    xsize(2) = 0.82*mri.dim(2)/xdim;
    ysize(1) = 0.82*mri.dim(3)/ydim;
    ysize(2) = 0.82*mri.dim(2)/ydim;
    
    %% create figure handles
    
    % axis handles will hold the anatomical data if present, along with labels etc.
    h1 = axes('position',[0.07 0.07+ysize(2)+0.05 xsize(1) ysize(1)]);
    h2 = axes('position',[0.07+xsize(1)+0.05 0.07+ysize(2)+0.05 xsize(2) ysize(1)]);
    h3 = axes('position',[0.07 0.07 xsize(1) ysize(2)]);
    set(h1,'Tag','ik','Visible','on','XAxisLocation','top');
    set(h2,'Tag','jk','Visible','on','XAxisLocation','top');
    set(h3,'Tag','ij','Visible','on');
    
    if hasIsotropicVoxels
      set(h1,'DataAspectRatio',[1 1 1]);
      set(h2,'DataAspectRatio',[1 1 1]);
      set(h3,'DataAspectRatio',[1 1 1]);
    end
    
    dat = double(mri.(cfg.parameter));
    dmin = min(dat(:));
    dmax = max(dat(:));
    dat  = (dat-dmin)./(dmax-dmin);
    
    x = 1:mri.dim(1);
    y = 1:mri.dim(2);
    z = 1:mri.dim(3);
    xc = round(mri.dim(1)/2);
    yc = round(mri.dim(2)/2);
    zc = round(mri.dim(3)/2);
    
    if isfield(cfg, 'pnt')
      pnt = cfg.pnt;
    else
      pnt = zeros(0,3);
    end
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
      '   a. select the position by clicking on it in any slice with the left mouse button\n',...
      '   b. identify it by pressing the letter corresponding to the fiducial/landmark:\n', fidexplanation1, fidexplanation2, ...
      '   You can mark the fiducials multiple times, until you are satisfied with the positions.\n',...
      '3. To change the display:\n',...
      '   a. press c on keyboard to toggle crosshair visibility\n',...
      '   b. press f on keyboard to toggle fiducial visibility\n',...
      '   c. press + or - on (numeric) keyboard to change the color range''s upper limit\n',...
      '4. To finalize markers and quit interactive mode, press q on keyboard\n'));
    
    % create structure to be passed to gui
    opt               = [];
    opt.dim           = mri.dim;
    opt.ijk           = [xc yc zc];
    opt.xsize         = xsize;
    opt.ysize         = ysize;
    opt.handlesaxes   = [h1 h2 h3];
    opt.handlesfigure = h;
    opt.quit          = false;
    opt.ana           = dat;
    opt.update        = [1 1 1];
    opt.init          = true;
    opt.tag           = 'ik';
    opt.mri           = mri;
    opt.showcrosshair = true;
    opt.showmarkers   = false;
    opt.markers       = {markerpos markerlabel markercolor};
    opt.clim          = [];
    opt.fiducial      = cfg.fiducial;
    opt.fidlabel      = fidlabel;
    opt.fidletter     = fidletter;
    opt.pnt           = pnt;
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
    
    % store the interactively determined fiducials in the configuration
    % the actual coordinate transformation will be done further down
    cfg.fiducial = opt.fiducial;
    
  case 'headshape'
    if strcmp(class(cfg.headshape), 'config')
      cfg.headshape = struct(cfg.headshape);
    end
    
    if ischar(cfg.headshape)
      % old-style specification, convert cfg into new representation
      cfg.headshape.headshape   = cfg.headshape;
      
      if isfield(cfg, 'scalpsmooth'),
        cfg.headshape.scalpsmooth = cfg.scalpsmooth;
        cfg = rmfield(cfg, 'scalpsmooth');
      end
      if isfield(cfg, 'scalpthreshold'),
        cfg.headshape.scalpthreshold = cfg.scalpthreshold;
        cfg = rmfield(cfg, 'scalpthreshold');
      end
    elseif isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
      % old-style specification, convert into new representation
      cfg.headshape.headshape   = cfg.headshape;
      if isfield(cfg, 'scalpsmooth'),
        cfg.headshape.scalpsmooth = cfg.scalpsmooth;
        cfg = rmfield(cfg, 'scalpsmooth');
      end
      if isfield(cfg, 'scalpthreshold'),
        cfg.headshape.scalpthreshold = cfg.scalpthreshold;
        cfg = rmfield(cfg, 'scalpthreshold');
      end
      
    elseif isstruct(cfg.headshape)
      % new-style specification, do nothing
    else
      error('incorrect specification of cfg.headshape');
    end
    if ischar(cfg.headshape.headshape)
      shape = ft_read_headshape(cfg.headshape.headshape);
    else
      shape = cfg.headshape.headshape;
    end
    shape = ft_convert_units(shape, 'mm');
    
    cfg.headshape.interactive    = ft_getopt(cfg.headshape, 'interactive', true);
    cfg.headshape.icp            = ft_getopt(cfg.headshape, 'icp',         true);
    cfg.headshape.scalpsmooth    = ft_getopt(cfg.headshape, 'scalpsmooth',    2, 1); % empty is OK
    cfg.headshape.scalpthreshold = ft_getopt(cfg.headshape, 'scalpthreshold', 0.1);
    
    dointeractive = istrue(cfg.headshape.interactive);
    doicp         = istrue(cfg.headshape.icp);
    
    if ~isfield(mri, 'scalp') || ~islogical(mri.scalp)
      % extract the scalp surface from the anatomical image
      tmpcfg        = [];
      tmpcfg.output = 'scalp';
      tmpcfg.scalpsmooth    = cfg.headshape.scalpsmooth;
      tmpcfg.scalpthreshold = cfg.headshape.scalpthreshold;
      if isfield(cfg, 'template')
        tmpcfg.template = cfg.template;
      end
      seg           = ft_volumesegment(tmpcfg, mri);
    else
      % use the scalp segmentation that is provided
      seg = mri;
    end
    
    tmpcfg             = [];
    tmpcfg.tissue      = 'scalp';
    tmpcfg.method      = 'projectmesh';%'isosurface';
    tmpcfg.numvertices = 20000;
    scalp              = ft_prepare_mesh(tmpcfg, seg);
    scalp              = ft_convert_units(scalp, 'mm');
    
    if dointeractive,
      tmpcfg                       = [];
      tmpcfg.template.elec         = shape;     % this is the Polhemus recorded headshape
      tmpcfg.template.elec.chanpos = shape.pnt;
      tmpcfg.individual.headshape  = scalp;     % this is the headshape extracted from the anatomical MRI
      tmpcfg.individual.headshapestyle = 'surface';
      tmpcfg = ft_interactiverealign(tmpcfg);
      M      = tmpcfg.m;
      cfg.transform_interactive = M;
      
      % touch it to survive trackconfig
      cfg.transform_interactive;
      
      % update the relevant geometrical info
      scalp  = ft_transform_geometry(M, scalp);
    end % dointeractive
    
    % always perform an icp-step, because this will give an estimate of the
    % initial distance of the corresponding points. depending on the value
    % for doicp, deal with the output differently
    if doicp,
      numiter = 50;
    else
      numiter = 1;
    end
    
    if ~isfield(cfg, 'weights')
      w = ones(size(shape.pnt,1),1);
    else
      w = cfg.weights(:);
      if numel(w)~=size(shape.pnt,1),
        error('number of weights should be equal to the number of points in the headshape');
      end
    end
    
    % the icp function wants this as a function handle.
    weights = @(x)assignweights(x,w);
    
    % construct the coregistration matrix
    nrm = normals(scalp.pnt, scalp.tri, 'vertex');
    [R, t, err, dummy, info] = icp(scalp.pnt', shape.pnt', numiter, 'Minimize', 'plane', 'Normals', nrm', 'Weight', weights, 'Extrapolation', true, 'WorstRejection', 0.05);
    
    if doicp,
      % create the additional transformation matrix and compute the
      % distance between the corresponding points, both prior and after icp
      
      % this one transforms from scalp 'headspace' to shape 'headspace'
      M2 = inv([R t;0 0 0 1]);
      
      % warp the extracted scalp points to the new positions
      scalp.pnt = ft_warp_apply(M2, scalp.pnt);
      
      target        = scalp;
      target.pos    = target.pnt;
      target.inside = (1:size(target.pos,1))';
      
      functional     = rmfield(shape,'pnt');
      functional.pow = info.distanceout(:);
      functional.pos = info.qout';
      
      tmpcfg              = [];
      tmpcfg.parameter    = 'pow';
      tmpcfg.interpmethod = 'sphere_avg';
      tmpcfg.sphereradius = 10;
      smoothdist          = ft_sourceinterpolate(tmpcfg, functional, target);
      scalp.distance      = smoothdist.pow(:);
      
      functional.pow      = info.distancein(:);
      smoothdist          = ft_sourceinterpolate(tmpcfg, functional, target);
      scalp.distancein    = smoothdist.pow(:);
      
      cfg.icpinfo = info;
      cfg.transform_icp = M2;
      
      % touch it to survive trackconfig
      cfg.icpinfo;
      cfg.transform_icp;
    else
      % compute the distance between the corresponding points, prior to
      % icp: this corresponds to the final result after interactive only
      
      M2 = eye(4); % this is needed later on
      
      target        = scalp;
      target.pos    = target.pnt;
      target.inside = (1:size(target.pos,1))';
      
      functional     = rmfield(shape,'pnt');
      functional.pow = info.distancein(:);
      functional.pos = info.qout';
      
      tmpcfg              = [];
      tmpcfg.parameter    = 'pow';
      tmpcfg.interpmethod = 'sphere_avg';
      tmpcfg.sphereradius = 10;
      smoothdist          = ft_sourceinterpolate(tmpcfg, functional, target);
      scalp.distance      = smoothdist.pow(:);
      
    end % doicp
    
    % create headshape structure for mri-based surface point cloud
    if isfield(mri, 'coordsys')
      scalp.coordsys = mri.coordsys;
      
      % coordsys is the same as input mri
      coordsys = mri.coordsys;
    else
      coordsys  = 'unknown';
    end
    
    % update the cfg
    cfg.headshape.headshape    = shape;
    cfg.headshape.headshapemri = scalp;
    
    % touch it to survive trackconfig
    cfg.headshape;
    
    if doicp && dointeractive
      transform = M2*M;
    elseif doicp
      transform = M2;
    elseif dointeractive
      transform = M;
    end
    
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
      
      vox = fopen(tmpname4);
      tmp = textscan(vox, '%f');
      fclose(vox);
      
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
    error('unsupported method "%s"', cfg.method);
end

if strcmp(cfg.method, 'fiducial') || strcmp(cfg.method, 'interactive')
  
  % the fiducial locations are specified in voxels, convert them to head
  % coordinates according to the existing transform matrix
  fid1_vox  = cfg.fiducial.(fidlabel{1});
  fid2_vox  = cfg.fiducial.(fidlabel{2});
  fid3_vox  = cfg.fiducial.(fidlabel{3});
  fid1_head = ft_warp_apply(mri.transform, fid1_vox);
  fid2_head = ft_warp_apply(mri.transform, fid2_vox);
  fid3_head = ft_warp_apply(mri.transform, fid3_vox);
  
  if length(fidlabel)>3
    % the 4th point is optional
    fid4_vox  = cfg.fiducial.(fidlabel{4});
    fid4_head = ft_warp_apply(mri.transform, fid4_vox);
  else
    fid4_head = [nan nan nan];
  end
  
  if ~any(isnan(fid4_head))
    [transform, coordsys] = ft_headcoordinates(fid1_head, fid2_head, fid3_head, fid4_head, cfg.coordsys);
  else
    [transform, coordsys] = ft_headcoordinates(fid1_head, fid2_head, fid3_head, cfg.coordsys);
  end
end

% copy the input anatomical or functional volume
realign = mri;

if ~isempty(transform) && ~any(isnan(transform(:)))
  % combine the additional transformation with the original one
  realign.transformorig = mri.transform;
  realign.transform     = transform * mri.transform;
  realign.coordsys      = coordsys;
else
  warning('no coordinate system realignment has been done');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous mri
ft_postamble history realign
ft_postamble savevar realign

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = assignweights(x, w)

% x is an indexing vector with the same number of arguments as w
y = w(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)
h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h,       'currentaxes');
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
  ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', opt.ijk, 'style', 'subplot', 'parents', [h1 h2 h3].*opt.update, 'doscale', false);
  
  opt.anahandles = findobj(opt.handlesfigure, 'type', 'surface')';
  parenttag  = get(cell2mat(get(opt.anahandles,'parent')),'tag');
  [i1,i2,i3] = intersect(parenttag, {'ik';'jk';'ij'});
  opt.anahandles = opt.anahandles(i3(i2)); % seems like swapping the order
  opt.anahandles = opt.anahandles(:)';
  set(opt.anahandles, 'tag', 'ana');
else
  ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', opt.ijk, 'style', 'subplot', 'surfhandle', opt.anahandles.*opt.update, 'doscale', false);
  
  if all(round([xi yi zi])<=mri.dim) && all(round([xi yi zi])>0)
    fprintf('==================================================================================\n');
    str = sprintf('voxel %d, index [%d %d %d]', sub2ind(mri.dim(1:3), round(xi), round(yi), round(zi)), round([xi yi zi]));
    
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
end

set(opt.handlesaxes(1),'Visible','on');
set(opt.handlesaxes(2),'Visible','on');
set(opt.handlesaxes(3),'Visible','on');

% make the last current axes current again
sel = findobj('type','axes','tag',tag);
if ~isempty(sel)
  set(opt.handlesfigure, 'currentaxes', sel(1));
end

if opt.init
  % draw the crosshairs for the first time
  hch1 = crosshair([xi 1 zi], 'parent', h1, 'color', 'yellow');
  hch3 = crosshair([xi yi opt.dim(3)], 'parent', h3, 'color', 'yellow');
  hch2 = crosshair([opt.dim(1) yi zi], 'parent', h2, 'color', 'yellow');
  opt.handlescross  = [hch1(:)';hch2(:)';hch3(:)'];
  opt.handlesmarker = [];
else
  % update the existing crosshairs, don't change the handles
  crosshair([xi 1 zi], 'handle', opt.handlescross(1, :));
  crosshair([opt.dim(1) yi zi], 'handle', opt.handlescross(2, :));
  crosshair([xi yi opt.dim(3)], 'handle', opt.handlescross(3, :));
end

if opt.showcrosshair
  set(opt.handlescross,'Visible','on');
else
  set(opt.handlescross,'Visible','off');
end

markercolor = {'r', 'g', 'b', 'y'};

delete(opt.handlesmarker(opt.handlesmarker(:)>0));
opt.handlesmarker = [];

for i=1:length(opt.fidlabel)
  pos = opt.fiducial.(opt.fidlabel{i});
  if any(isnan(pos))
    continue
  end
  
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

if opt.showmarkers
  set(opt.handlesmarker,'Visible','on');
else
  set(opt.handlesmarker,'Visible','off');
end

if opt.init
  opt.init = false;
end

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

curr_ax = get(h,       'currentaxes');
tag     = get(curr_ax, 'tag');

if isempty(key)
  % this happens if you press the apple key
  key = '';
end

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
    sel = strcmp(key, opt.fidletter);
    fprintf('==================================================================================\n');
    fprintf('selected %s\n', opt.fidlabel{sel});
    opt.fiducial.(opt.fidlabel{sel}) = opt.ijk;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
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
    opt.clim(2) = opt.clim(2)-cscalefactor;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
  case {45 'shift+hyphen'} % numpad -
    if isempty(opt.clim)
      opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
    end
    % increase color scale range by 10%
    cscalefactor = (opt.clim(2)-opt.clim(1))/10;
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

uiresume(h);

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

uiresume;

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
uiresume;

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

setappdata(h, 'opt', opt);
uiresume;

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
    key = eventdata.Character;
  end
end

if ~isempty(eventdata.Modifier)
  key = [eventdata.Modifier{1} '+' key];
end
