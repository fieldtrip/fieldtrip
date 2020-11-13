function [sourcemodel, cfg] = ft_prepare_sourcemodel(cfg)

% FT_PREPARE_SOURCEMODEL constructs a source model, for example a 3D grid or a
% cortical sheet. The source model that can be used for source reconstruction,
% beamformer scanning, linear estimation and MEG interpolation.
%
% Use as
%   sourcemodel = ft_prepare_sourcemodel(cfg)
% where the details of the configuration structure determine how the source
% model will be constructed.
%
% The different approaches for constructing a source model are
%   cfg.method = 'basedongrid'        regular 3D grid with explicit specification
%                'basedonpos'         regular 3D grid with specification of the resolution
%                'basedonshape'       surface mesh based on inward shifted head surface from external file
%                'basedonmri'         regular 3D grid, based on segmented MRI, restricted to gray matter
%                'basedonmni'         regular 3D grid, based on a warped template grid, based on the MNI brain
%                'basedoncortex'      cortical sheet from external software such as Caret or FreeSurfer, can also be two separate hemispheres
%                'basedonresolution'  regular 3D grid with specification of the resolution
%                'basedonvol'         surface mesh based on inward shifted brain surface from volume conductor
%                'basedonfile'        the sourcemodel should be read from file
%                'basedoncentroids'   irregular 3D grid based on volumetric mesh
% The default for cfg.method is to determine the approach automatically, based on
% the configuration options that you specify.
%
% BASEDONRESOLUTION - uses an explicitly specified grid, or with the desired
% resolution, according to the following configuration options:
%   cfg.xgrid         = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.ygrid         = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.zgrid         = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.resolution    = number (e.g. 1 cm) for automatic grid generation
%
% BASEDONPOS - places sources on positions that you explicitly specify,
% according to the following configuration options:
%   cfg.sourcemodel.pos       = N*3 matrix with position of each source
%   cfg.sourcemodel.inside    = N*1 vector with boolean value whether position is inside brain (optional)
%   cfg.sourcemodel.dim       = [Nx Ny Nz] vector with dimensions in case of 3D grid (optional)
% The following fields (from FT_PRERARE_LEADFIELD or FT_SOURCEANALYSIS) are
% not used in this function, but will be copied along to the output:
%   cfg.sourcemodel.leadfield = cell-array
%   cfg.sourcemodel.filter    = cell-array
%   cfg.sourcemodel.subspace
%   cfg.sourcemodel.lbex
%
% BASEDONMNI - uses source positions from a template sourcemodel that is
% inversely warped from MNI coordinates to the individual subjects MRI.
% It uses the following configuration options:
%   cfg.mri             = structure with anatomical MRI model or filename, see FT_READ_MRI
%   cfg.nonlinear       = 'no' (or 'yes'), use non-linear normalization
%   cfg.resolution      = number (e.g. 6) of the resolution of the template MNI grid, defined in mm
%   cfg.template        = specification of a template sourcemodel as structure, or the filename of a template sourcemodel (defined in MNI space)
% Either cfg.resolution or cfg.template needs to be defined; if both are defined, cfg.template prevails.
%
% BASEDONMRI - makes a segmentation of the individual anatomical MRI and places
% sources in the grey matter. It uses the following configuration options:
%   cfg.mri             = can be filename, MRI structure or segmented MRI structure
%   cfg.threshold       = 0.1, relative to the maximum value in the segmentation
%   cfg.smooth          = 5, smoothing in voxels
%
% BASEDONCORTEX - places sources on the vertices of a cortical surface description
%   cfg.headshape       = string, should be a *.fif file
%
% BASEDONCENTROIDS - places sources on the centroids of a volumetric mesh
%   cfg.headmodel       = volumetric mesh
%   cfg.headmodel.type  = 'simbio';
%
% Other configuration options include
%   cfg.unit            = string, can be 'mm', 'cm', 'm' (default is automatic)
%   cfg.tight           = 'yes' or 'no' (default is automatic)
%   cfg.inwardshift     = number, amount to shift the innermost surface of the headmodel inward when determining
%                         whether sources are inside or outside the source compartment (default = 0)
%   cfg.moveinward      = number, amount to move sources inward to ensure a certain minimal distance to the innermost
%                         surface of the headmodel (default = 0)
%   cfg.movetocentroids = 'yes' or 'no', move the dipoles to the centroids of the hexahedral
%                         or tetrahedral mesh (default = 'no')
%   cfg.spherify        = 'yes' or 'no', scale the source model so that it fits inside a sperical
%                         volume conduction model (default = 'no')
%   cfg.symmetry        = 'x', 'y' or 'z' symmetry for two dipoles, can be empty (default = [])
%   cfg.headshape       = a filename for the headshape, a structure containing a single surface,
%                         or a Nx3 matrix with headshape surface points (default = [])
%   cfg.spmversion      = string, 'spm2', 'spm8', 'spm12' (default = 'spm12')
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad          = structure with gradiometer definition or filename, see FT_READ_SENS
%
% The headmodel or volume conduction model can be specified as
%   cfg.headmodel     = structure with volume conduction model or filename, see FT_PREPARE_HEADMODEL
%
% The cfg.inwardshift option can be used for 3D grids to specify a positive (inward)
% or negative (outward) number to shift the innermost surface of the headmodel
% (usually the skull) when determining whether sources are to be flagged as inside or
% outside the source compartment. Only sources flagged as inside will be considered
% for subsequent source reconstructions. An ourward shift can be useful for a
% spherical or singleshell MEG headmodel. For a source model based on a cortical
% sheet in general you want all sources to be considered inside. For a BEM headmodel
% (EEG or MEG), there should never be any sources outside the actual source
% compartment.
%
% The cfg.moveinward option can be used for a source model based on a cortical sheet
% to push the sources inward a little bit to ensure sufficient distance to the
% innermost surface of a BEM headmodel (EEG or MEG).
%
% See also FT_PREPARE_LEADFIELD, FT_PREPARE_HEADMODEL, FT_SOURCEANALYSIS,
% FT_DIPOLEFITTING, FT_MEGREALIGN

% Copyright (C) 2004-2020, Robert Oostenveld
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
ft_preamble provenance headmodel sens
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', 'mriunits');
cfg = ft_checkconfig(cfg, 'forbidden', 'sourceunits');
cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'grid',    'sourcemodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed', {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed', {'optofile', 'opto'});
cfg = ft_checkconfig(cfg, 'renamedval', {'unit', 'auto', []});

cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'});  % this is moved to cfg.sourcemodel.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.unit by the subsequent createsubcfg

% put the low-level options pertaining to the sourcemodel in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg', {'sourcemodel'});
% move some fields from cfg.sourcemodel back to the top-level configuration
cfg = ft_checkconfig(cfg, 'createtopcfg', {'sourcemodel'});

% set the defaults
cfg.moveinward        = ft_getopt(cfg, 'moveinward'); % the default is automatic and depends on a triangulation being present
cfg.spherify          = ft_getopt(cfg, 'spherify', 'no');
cfg.headshape         = ft_getopt(cfg, 'headshape');
cfg.symmetry          = ft_getopt(cfg, 'symmetry');
cfg.spmversion        = ft_getopt(cfg, 'spmversion', 'spm12');
cfg.headmodel         = ft_getopt(cfg, 'headmodel');
cfg.sourcemodel       = ft_getopt(cfg, 'sourcemodel');
cfg.unit              = ft_getopt(cfg, 'unit');
cfg.method            = ft_getopt(cfg, 'method'); % the default is to do automatic detection further down
cfg.movetocentroids   = ft_getopt(cfg, 'movetocentroids', 'no');

% this option was deprecated on 12 Aug 2020
if isfield(cfg, 'warpmni')
  % prior to the introduction of cfg.method we used cfg.warpmni to separate between
  % basedonmni and basedonmri, which both require cfg.mri to be present
  if isfield(cfg, 'mri') && istrue(cfg.warpmni)
    ft_warning('please specify cfg.method=''basedonmni'' instead of cfg.warpmni=''yes''');
    cfg.method = 'basedonmni';
  elseif isfield(cfg, 'mri') && ~istrue(cfg.warpmni)
    ft_warning('please specify cfg.method=''basedonmri'' instead of cfg.warpmni=''no''');
    cfg.method = 'basedonmri';
  end
  cfg = rmfield(cfg, 'warpmni');
end

% this code expects the inside to be represented as a logical array
cfg = ft_checkconfig(cfg, 'inside2logical', 'yes');
if isfield(cfg, 'sourcemodel')
  cfg.sourcemodel = ft_checkconfig(cfg.sourcemodel, 'renamed',  {'pnt' 'pos'});
end
if isfield(cfg, 'template')
  cfg.template = ft_checkconfig(cfg.template, 'renamed',  {'pnt' 'pos'});
end

if isfield(cfg, 'resolution') && isfield(cfg, 'xgrid') && ~ischar(cfg.xgrid)
  ft_error('You cannot specify cfg.resolution and an explicit cfg.xgrid simultaneously');
end
if isfield(cfg, 'resolution') && isfield(cfg, 'ygrid') && ~ischar(cfg.ygrid)
  ft_error('You cannot specify cfg.resolution and an explicit cfg.ygrid simultaneously');
end
if isfield(cfg, 'resolution') && isfield(cfg, 'zgrid') && ~ischar(cfg.zgrid)
  ft_error('You cannot specify cfg.resolution and an explicit cfg.zgrid simultaneously');
end

% the source model can be constructed in a number of ways
if isempty(cfg.method)
  if isfield(cfg, 'sourcemodel') && ischar(cfg.sourcemodel)
    cfg.method = 'basedonfile';
  elseif isfield(cfg, 'xgrid') && ~ischar(cfg.xgrid)
    cfg.method = 'basedongrid'; % regular 3D grid with explicit specification
  elseif isfield(cfg.sourcemodel, 'pos')
    cfg.method = 'basedonpos'; % using user-supplied positions, which can be regular or irregular
  elseif ~isempty(cfg.headshape)
    cfg.method = 'basedonshape'; % surface mesh based on inward shifted head surface from external file
  elseif isfield(cfg, 'mri')
    cfg.method = 'basedonmri'; % regular 3D grid, based on segmented MRI, restricted to gray matter
  elseif isfield(cfg, 'headshape') && (iscell(cfg.headshape) || any(ft_filetype(cfg.headshape, {'neuromag_fif', 'freesurfer_triangle_binary', 'caret_surf', 'gifti'})))
    cfg.method = 'basedoncortex'; % cortical sheet from external software such as Caret or FreeSurfer, can also be two separate hemispheres
  elseif isfield(cfg, 'resolution')
    cfg.method = 'basedonresolution'; % regular 3D grid with specification of the resolution
  elseif ~isempty(cfg.headmodel)
    cfg.method = 'basedonvol'; % surface mesh based on inward shifted brain surface from volume conductor
  else
    ft_error('incorrect cfg specification for constructing a sourcemodel');
  end
else
  cfg = ft_checkconfig(cfg, 'allowedval', {'method', 'basedongrid', 'basedonpos', 'basedonshape', ...
    'basedonmri', 'basedonmni', 'basedoncortex', 'basedonresolution', 'basedonvol', 'basedonfile','basedoncentroids'});
end

% these are mutually exclusive, but printing all requested methods here
% facilitates debugging of weird configs. Also specify the defaults here to
% keep the overview
switch cfg.method
  case 'basedonfile'
    fprintf('reading sourcemodel from file\n');
    cfg.tight       = ft_getopt(cfg, 'tight',   'no');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection
    
  case 'basedonresolution'
    fprintf('creating sourcemodel based on automatic 3D grid with specified resolution\n');
    cfg.xgrid       = ft_getopt(cfg, 'xgrid',  'auto');
    cfg.ygrid       = ft_getopt(cfg, 'ygrid',  'auto');
    cfg.zgrid       = ft_getopt(cfg, 'zgrid',  'auto');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection
    cfg.tight       = ft_getopt(cfg, 'tight',   'yes');
    
  case 'basedongrid'
    fprintf('creating sourcemodel based on user specified 3D grid\n');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection
    cfg.tight       = ft_getopt(cfg, 'tight',   'yes');
    
  case 'basedonpos'
    fprintf('creating sourcemodel based on user specified dipole positions\n');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection
    cfg.tight       = ft_getopt(cfg, 'tight',    'no');
    
  case 'basedonshape'
    fprintf('creating sourcemodel based on inward-shifted head shape\n');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift',  0); % in this case for inside detection
    cfg.spheremesh  = ft_getopt(cfg, 'spheremesh', 642);
    cfg.tight       = ft_getopt(cfg, 'tight',    'yes');
    
  case 'basedoncortex'
    cfg.tight       = ft_getopt(cfg, 'tight', 'yes');
    
  case 'basedonmri'
    fprintf('creating sourcemodel based on an anatomical volume\n');
    cfg.threshold   = ft_getopt(cfg, 'threshold', 0.1); % relative
    cfg.smooth      = ft_getopt(cfg, 'smooth',      5); % in voxels
    cfg.tight       = ft_getopt(cfg, 'tight',   'yes');
    
  case 'basedonvol'
    fprintf('creating sourcemodel based on inward-shifted brain surface from volume conductor model\n');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift',   0); % in this case for inside detection
    cfg.spheremesh  = ft_getopt(cfg, 'spheremesh',  642);
    cfg.tight       = ft_getopt(cfg, 'tight',      'no');
    
  case 'basedonmni'
    cfg.tight       = ft_getopt(cfg.sourcemodel, 'tight',       'no');
    cfg.nonlinear   = ft_getopt(cfg.sourcemodel, 'nonlinear',   'no');
    
  case 'basedoncentroids'
    fprintf('creating sourcemodel based on volumetric mesh centroids\n');
    cfg.tight       = ft_getopt(cfg.sourcemodel, 'tight',       'no');
    cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 0); % in this case for inside detection
end

if (isfield(cfg, 'smooth') && ~strcmp(cfg.smooth, 'no')) || strcmp(cfg.method, 'basedonmni')
  % check that the preferred SPM version is on the path
  ft_hastoolbox(cfg.spmversion, 1);
end

% start with an empty structure
sourcemodel = [];

% get the volume conduction model
if ischar(cfg.headmodel)
  headmodel = ft_read_headmodel(cfg.headmodel);
else
  % ensure that the volume conduction model is up-to-date
  headmodel = ft_datatype_headmodel(cfg.headmodel);
end

% get the gradiometer or electrode definition
try
  sens = ft_fetch_sens(cfg);
catch
  sens = [];
end

if strcmp(cfg.method, 'basedonfile')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the source model from a MATLAB file
  % this needs to be done prior to determining the default units
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cfg.sourcemodel = loadvar(cfg.sourcemodel, 'sourcemodel');
  sourcemodel = cfg.sourcemodel;
end

if isempty(cfg.unit)
  if isfield(cfg.sourcemodel, 'unit') && ~isempty(cfg.sourcemodel.unit)
    % take the existing source model units
    cfg.unit = cfg.sourcemodel.unit;
  elseif isfield(cfg.sourcemodel, 'pos') && size(cfg.sourcemodel.pos,1)>10
    % estimate the units based on the existing source positions
    cfg.sourcemodel = ft_determine_units(cfg.sourcemodel);
    cfg.unit = cfg.sourcemodel.unit;
  elseif strcmp(cfg.method, 'basedonmni') && ~isempty(cfg.mri.unit)
    % take the existing MRI units
    cfg.unit = cfg.mri.unit;
  elseif ~isempty(sens)
    % take the units from the gradiometer or electrode array
    cfg.unit = sens.unit;
  elseif ~isempty(headmodel)
    % take the units from the volume conduction model
    cfg.unit = headmodel.unit;
  else
    ft_warning('assuming "cm" as default units for source model');
    cfg.unit = 'cm';
  end
end

% convert the sensor array to the desired units for the source model
if ~isempty(sens)
  sens = ft_convert_units(sens, cfg.unit);
end

% convert the head model to the desired units for the source model
if ~isempty(headmodel)
  headmodel = ft_convert_units(headmodel, cfg.unit);
end

switch cfg.method
  case 'basedonresolution'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct a regular 3D grid that spans a box encompassing all electrode
    % or gradiometer coils, this will typically also cover the complete brain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(sens) && isfield(headmodel, 'chanpos')
      % determine the bounding box of the sensor array
      minpos = min(sens.chanpos,[],1);
      maxpos = max(sens.chanpos,[],1);
    elseif ~isempty(headmodel)
      % determine the bounding box of the volume conduction model
      if isfield(headmodel, 'bnd') && isfield(headmodel.bnd, 'pos')
        pos = cat(1, headmodel.bnd(:).pos);
      elseif isfield(headmodel, 'bnd') && isfield(headmodel.bnd, 'pnt')
        pos = cat(1, headmodel.bnd(:).pnt);
      elseif isfield(headmodel, 'pos')
        pos = headmodel.pos;
      elseif ft_headmodeltype(headmodel, 'localspheres')
        pos = headsurface(headmodel, sens);
      elseif ft_headmodeltype(headmodel, 'singlesphere')
        pos = [
          headmodel.o - headmodel.r
          headmodel.o + headmodel.r
          ];
      elseif ft_headmodeltype(headmodel, 'concentricspheres')
        pos = [
          headmodel.o - max(headmodel.r)
          headmodel.o + max(headmodel.r)
          ];
      end
      minpos = min(pos,[],1);
      maxpos = max(pos,[],1);
      % add a few percent on either side
      minpos(minpos<0) = minpos(minpos<0).*1.08;
      maxpos(maxpos>0) = maxpos(maxpos>0).*1.08;
      minpos(minpos>0) = minpos(minpos>0).*0.92;
      maxpos(maxpos<0) = maxpos(maxpos<0).*0.92;
    else
      ft_error('creating a 3D grid based on resolution requires either sensor positions or a headmodel to estimate the extent');
    end
    
    fprintf('creating 3D grid with %g %s resolution\n', cfg.resolution, cfg.unit);
    
    % round the bounding box limits to the nearest cm
    switch cfg.unit
      case 'm'
        minpos = floor(minpos*100)/100;
        maxpos = ceil(maxpos*100)/100;
      case 'cm'
        minpos = floor(minpos);
        maxpos = ceil(maxpos);
      case 'mm'
        minpos = floor(minpos/10)*10;
        maxpos = ceil(maxpos/10)*10;
    end
    
    if ischar(cfg.xgrid) && strcmp(cfg.xgrid, 'auto')
      xgrid = minpos(1):cfg.resolution:maxpos(1);
    end
    if ischar(cfg.ygrid) && strcmp(cfg.ygrid, 'auto')
      ygrid = minpos(2):cfg.resolution:maxpos(2);
    end
    if ischar(cfg.zgrid) && strcmp(cfg.zgrid, 'auto')
      zgrid = minpos(3):cfg.resolution:maxpos(3);
    end
    sourcemodel.dim   = [length(xgrid) length(ygrid) length(zgrid)];
    [X, Y, Z]  = ndgrid(xgrid, ygrid, zgrid);
    sourcemodel.pos   = [X(:) Y(:) Z(:)];
    sourcemodel.unit  = cfg.unit;
    fprintf('initial 3D grid dimensions are [%d %d %d]\n', sourcemodel.dim(1), sourcemodel.dim(2), sourcemodel.dim(3));
    
  case 'basedongrid'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % a detailed xgrid/ygrid/zgrid has been specified, the other details
    % still need to be determined
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sourcemodel.dim   = [length(cfg.xgrid) length(cfg.ygrid) length(cfg.zgrid)];
    [X, Y, Z]  = ndgrid(cfg.xgrid, cfg.ygrid, cfg.zgrid);
    sourcemodel.pos   = [X(:) Y(:) Z(:)];
    sourcemodel.unit  = cfg.unit;
    
  case 'basedonpos'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % source positions are already specified in the configuration, reuse as much of the
    % prespecified model as possible (but only known objects)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sourcemodel = keepfields(cfg.sourcemodel, {'pos', 'unit', 'xgrid', 'ygrid', 'zgrid', 'mom', 'tri', 'dim', 'transform', 'inside', 'lbex', 'subspace', 'leadfield', 'filter', 'label', 'leadfielddimord'});
    
  case 'basedonmri'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct a grid based on the segmented MRI that is provided in the
    % configuration, only voxels in gray matter will be used
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ischar(cfg.mri)
      mri = ft_read_mri(cfg.mri);
    else
      mri = cfg.mri;
    end
    
    % ensure the mri to have units
    if ~isfield(mri, 'unit')
      mri = ft_determine_units(mri);
    end
    
    if ~isfield(cfg, 'resolution')
      switch cfg.unit
        case 'mm'
          cfg.resolution = 10;
        case 'cm'
          cfg.resolution = 1;
        case 'dm'
          cfg.resolution = 0.1;
        case 'm'
          cfg.resolution = 0.01;
      end
    end
    
    issegmentation = false;
    if isfield(mri, 'gray')
      % this is not a boolean segmentation, but based on tissue probability
      % maps, being the original implementation here.
      dat = double(mri.gray);
      
      % apply a smoothing of a certain amount of voxels
      if ~strcmp(cfg.smooth, 'no')
        dat = volumesmooth(dat, cfg.smooth, 'MRI gray matter');
      end
      
    elseif isfield(mri, 'anatomy')
      % this could be a tpm stored on disk, i.e. the result of
      % ft_volumesegment. Reading it in always leads to the field 'anatomy'.
      % Note this could be any anatomical mask
      dat = double(mri.anatomy);
      
      % apply a smoothing of a certain amount of voxels
      if ~strcmp(cfg.smooth, 'no')
        dat = volumesmooth(dat, cfg.smooth, 'anatomy');
      end
      
    elseif ft_datatype(mri, 'segmentation')
      % this is a proper segmentation, where a set of boolean masks is in the
      % input, or and indexed volume, along with labels. FIXME for now still
      % only works for boolean volumes.
      issegmentation = true;
      fn = booleanfields(mri);
      if isempty(fn)
        % convert indexed segmentation into probabilistic
        mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic');
        fn  = booleanfields(mri);
      end
      
      dat = false(mri.dim);
      for i=1:numel(fn)
        if ~strcmp(cfg.smooth, 'no')
          mri.(fn{i}) = volumesmooth(double(mri.(fn{i})), cfg.smooth, fn{i}) > cfg.threshold;
        end
        dat = dat | mri.(fn{i});
      end
      dat = double(dat);
    else
      ft_error('cannot determine the format of the segmentation in cfg.mri');
    end
    
    % determine for each voxel whether it belongs to the grey matter
    fprintf('thresholding MRI data at a relative value of %f\n', cfg.threshold);
    head = dat./max(dat(:)) > cfg.threshold;
    
    % convert the source/functional data into the same units as the anatomical MRI
    scale = ft_scalingfactor(cfg.unit, mri.unit);
    
    ind                 = find(head(:));
    fprintf('%d from %d voxels in the segmentation are marked as ''inside'' (%.0f%%)\n', length(ind), numel(head), 100*length(ind)/numel(head));
    [X,Y,Z]             = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));             % create the grid in MRI-coordinates
    posmri              = [X(ind) Y(ind) Z(ind)];                                       % take only the inside voxels
    poshead             = ft_warp_apply(mri.transform, posmri);                         % transform to head coordinates
    resolution          = cfg.resolution*scale;                                         % source and mri can be expressed in different units (e.g. cm and mm)
    xgrid               = floor(min(poshead(:,1))):resolution:ceil(max(poshead(:,1)));  % create the grid in head-coordinates
    ygrid               = floor(min(poshead(:,2))):resolution:ceil(max(poshead(:,2)));  % with 'consistent' x,y,z definitions
    zgrid               = floor(min(poshead(:,3))):resolution:ceil(max(poshead(:,3)));
    [X,Y,Z]             = ndgrid(xgrid,ygrid,zgrid);
    pos2head            = [X(:) Y(:) Z(:)];
    pos2mri             = ft_warp_apply(inv(mri.transform), pos2head);                  % transform to MRI voxel coordinates
    pos2mri             = round(pos2mri);
    inside              = getinside(pos2mri, head);                                     % use helper subfunction
    
    sourcemodel.pos     = pos2head/scale;                                               % convert to source units
    sourcemodel.dim     = [length(xgrid) length(ygrid) length(zgrid)];
    sourcemodel.inside  = inside(:);
    sourcemodel.unit    = cfg.unit;
    
    if issegmentation
      % pass on the segmentation information on the grid points, the
      % individual masks have been smoothed above
      fn = booleanfields(mri);
      for i=1:numel(fn)
        sourcemodel.(fn{i}) = getinside(pos2mri, mri.(fn{i}));
      end
      % convert back is not in general possible because the masks can be
      % overlapping due to smoothing
      % sourcemodel = ft_datatype_segmentation(sourcemodel, 'segmentationstyle', segstyle);
    end
    fprintf('the full grid contains %d points\n', numel(sourcemodel.inside));
    fprintf('%d grid points are marked as inside the brain\n',  sum(sourcemodel.inside));
    
  case 'basedoncortex'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read it from a *.fif file that was created using Freesurfer and MNE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iscell(cfg.headshape)
      % FIXME loop over all files, this should be two hemispheres
      keyboard
    else
      shape = ft_read_headshape(cfg.headshape);
    end
    % ensure that the headshape is in the same units as the source
    shape     = ft_convert_units(shape, cfg.unit);
    % return both the vertices and triangles from the cortical sheet
    sourcemodel.pos  = shape.pos;
    sourcemodel.tri  = shape.tri;
    sourcemodel.unit = shape.unit;
    
  case 'basedonshape'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use the headshape  to make a superficial dipole layer (e.g.
    % for megrealign). Assume that all points are inside the volume.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the surface describing the head shape
    [headshape.pos, headshape.tri] = headsurface([], [], 'headshape', cfg.headshape);
    
    % ensure that the headshape is in the same units as the source model
    headshape = ft_convert_units(headshape, cfg.unit);
    
    % note that cfg.inwardshift should be expressed in the units consistent with the data
    sourcemodel.pos     = headsurface([], [], 'headshape', headshape, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
    sourcemodel.tri     = headshape.tri;
    sourcemodel.unit    = headshape.unit;
    sourcemodel.inside  = true(size(sourcemodel.pos,1),1);
    
  case 'basedonvol'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use the volume conduction model to make a superficial dipole layer (e.g.
    % for megrealign). Assume that all points are inside the volume.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % please note that cfg.inwardshift should be expressed in the units consistent with cfg.unit
    sourcemodel.pos     = headsurface(headmodel, sens, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
    sourcemodel.unit    = cfg.unit;
    sourcemodel.inside  = true(size(sourcemodel.pos,1),1);
    
  case 'basedonmni'
    if ~isfield(cfg, 'template') && ~isfield(cfg, 'resolution')
      ft_error('you either need to specify the filename of a template grid in cfg.template, or a resolution in cfg.resolution');
    elseif isfield(cfg, 'template')
      % let the template filename prevail
      fname = cfg.template;
    elseif isfield(cfg, 'resolution')
      % use one of the templates that are in Fieldtrip, this requires a resolution
      if isequal(cfg.resolution, round(cfg.resolution))
        fname = sprintf('standard_sourcemodel3d%dmm.mat', cfg.resolution);
      else
        fname = sprintf('standard_sourcemodel3d%dpoint%dmm.mat', floor(cfg.resolution), round(10*(cfg.resolution-floor(cfg.resolution))));
      end
      if ~exist(fname, 'file')
        ft_error('the MNI template grid based on the specified resolution does not exist');
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check whether the mni template grid exists for the specified resolution
    % if not create it: FIXME (this needs to be done still)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get the mri
    if ischar(cfg.mri)
      mri = ft_read_mri(cfg.mri);
    else
      mri = cfg.mri;
    end
    
    % get the template grid
    if ischar(fname)
      mnigrid = loadvar(fname, 'sourcemodel');
    else
      mnigrid = cfg.template;
    end
    
    % ensure these to have units in mm, the conversion of the source model is done further down
    mri     = ft_convert_units(mri,     'mm');
    mnigrid = ft_convert_units(mnigrid, 'mm');
    
    % ensure that it is specified with logical inside
    mnigrid = fixinside(mnigrid);
    
    % spatial normalisation of mri and construction of subject specific sourcemodel positions
    tmpcfg           = keepfields(cfg, {'spmversion', 'spmmethod', 'nonlinear'});
    if isfield(cfg, 'templatemri')
      % this option is called differently for the two functions
      tmpcfg.template = cfg.templatemri;
    end
    normalise = ft_volumenormalise(tmpcfg, mri);
    
    if ~isfield(normalise, 'params') && ~isfield(normalise, 'initial')
      % this is for older implementations of FT_VOLUMENORMALISE
      fprintf('applying an inverse warp based on a linear transformation only\n');
      sourcemodel.pos = ft_warp_apply(inv(normalise.cfg.final), mnigrid.pos);
    else
      % the normalisation from original subject head coordinates to MNI consists of an initial linear, followed by a nonlinear transformation
      % the reverse transformations need to be done to get from MNI to the original subject head coordinates
      % first apply the inverse of the nonlinear transformation, followed by the inverse of the initial linear transformation
      sourcemodel.pos = ft_warp_apply(inv(normalise.initial), ft_warp_apply(normalise.params, mnigrid.pos, 'sn2individual'));
    end
    % copy some of the fields over from the input arguments
    sourcemodel = copyfields(mri,       sourcemodel, {'unit', 'coordsys'});
    sourcemodel = copyfields(mnigrid,   sourcemodel, {'dim', 'tri', 'inside'});
    sourcemodel = copyfields(normalise, sourcemodel, {'params', 'initial'});
    if ft_datatype(mnigrid, 'parcellation')
      % copy the boolean fields over from the template MNI grid
      sourcemodel = copyfields(mnigrid, sourcemodel, booleanfields(mnigrid));
    end
    
  case 'basedoncentroids'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the centroids of each volume element of a FEM mesh
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this will also copy some fields from the headmodel, such as tissue, tissuelabel, unit, coordsys
    sourcemodel = compute_centroids(headmodel);
end

if isfield(sourcemodel, 'unit')
  % in most cases the source model will already be in the desired units, but e.g. for
  % "basedonmni" it will be in 'mm' since in the spatial normalization the MRI and
  % the MNI template are converted to 'mm'
  sourcemodel = ft_convert_units(sourcemodel, cfg.unit);
else
  % the units were specified by the user or determined automatically, assign them to the source model
  sourcemodel.unit = cfg.unit;
end

% do some sanity checks
if isfield(sourcemodel, 'filter')
  assert(numel(sourcemodel.filter) == size(sourcemodel.pos, 1), 'the number of precomputed filters does not match number of source positions');
end
if isfield(sourcemodel, 'leadfield')
  assert(numel(sourcemodel.leadfield) == size(sourcemodel.pos, 1), 'the number of precomputed leadfields does not match number of source positions');
end

if strcmp(cfg.spherify, 'yes')
  if ~ft_headmodeltype(headmodel, 'singlesphere') && ~ft_headmodeltype(headmodel, 'concentricspheres')
    ft_error('this only works for spherical volume conduction models');
  end
  % deform the cortex so that it fits in a unit sphere
  pos = mesh_spherify(sourcemodel.pos, [], 'shift', 'range');
  % scale it to the radius of the innermost sphere, make it a tiny bit smaller to
  % ensure that the support point with the exact radius 1 is still inside the sphere
  pos = pos*min(headmodel.r)*0.999;
  pos(:,1) = pos(:,1) + headmodel.o(1);
  pos(:,2) = pos(:,2) + headmodel.o(2);
  pos(:,3) = pos(:,3) + headmodel.o(3);
  sourcemodel.pos = pos;
end

if ~isempty(cfg.moveinward)
  if ~ismember(cfg.method, {'basedonshape', 'basedoncortex', 'basedonvol', 'basedonfile'})
    ft_warning('cfg.moveinward is designed to work with surface based sourcemodels, not with 3D grid sourcemodels.')
  end
  % construct a triangulated boundary of the source compartment
  [pos1, tri1] = headsurface(headmodel, [], 'inwardshift', cfg.moveinward, 'surface', 'brain');
  inside = bounding_mesh(sourcemodel.pos, pos1, tri1);
  if ~all(inside)
    pos2 = sourcemodel.pos(~inside,:);
    [dum, pos3] = project_elec(pos2, pos1, tri1);
    sourcemodel.pos(~inside,:) = pos3;
  end
  if cfg.moveinward>cfg.inwardshift
    sourcemodel.inside  = true(size(sourcemodel.pos,1),1);
  end
end

if strcmp(cfg.movetocentroids, 'yes')
  % compute centroids of the tetrahedral or hexahedral mesh
  centroids = compute_centroids(headmodel);
  
  % move the dipole positions in the sourcemodel to the closest centroid
  grid_shifted = zeros(size(sourcemodel.pos));
  for i = 1:length(sourcemodel.pos)
    [dum, amin] = min(sum((sourcemodel.pos(i,:) - centroids.pos).^2,2));
    grid_shifted(i,:) = centroids.pos(amin,:);
  end
  % eliminate duplicates, this applies for example if cfg.resolution is smaller than the mesh resolution
  sourcemodel.pos = unique(grid_shifted, 'rows', 'stable');
  
  % the positions are not on a regular 3D grid any more, hence dim does not apply
  sourcemodel = removefields(sourcemodel, {'dim'});
end

% determine the dipole locations that are inside the source compartment of the
% volume conduction model, i.e. inside the brain
if ~isfield(sourcemodel, 'inside')
  sourcemodel.inside = ft_inside_headmodel(sourcemodel.pos, headmodel, 'grad', sens, 'headshape', cfg.headshape, 'inwardshift', cfg.inwardshift); % this returns a boolean vector
else
  if isfield(cfg, 'inwardshift') && isfield(cfg, 'template')
    % warn about inwardshift not having an effect as inside is already specified as well
    % warning should only be issued for templates, inwardshift can also be present for surface meshes
    ft_warning('Inside dipole locations already determined by a template, cfg.inwardshift has no effect.')
  end
end

if strcmp(cfg.tight, 'yes')
  if ~isfield(sourcemodel, 'dim')
    ft_error('tightgrid only works for positions on a regular 3D grid');
  end
  fprintf('%d dipoles inside, %d dipoles outside brain\n', sum(sourcemodel.inside), sum(~sourcemodel.inside));
  fprintf('making tight grid\n');
  boolvol = reshape(sourcemodel.inside, sourcemodel.dim);
  xsel    = squeeze(sum(sum(boolvol,3),2))>0;
  ysel    = squeeze(sum(sum(boolvol,3),1))>0;
  zsel    = squeeze(sum(sum(boolvol,2),1))>0;
  boolvol(xsel,ysel,zsel) = true; % update the volume to contain the to-be-selected entries
  sel     = boolvol(:);
  
  % update the boolean fields, this requires the original dim
  fn = booleanfields(sourcemodel);
  for i=1:numel(fn)
    sourcemodel.(fn{i}) = sourcemodel.(fn{i})(sel);
  end
  
  % update the grid locations that are marked as inside the brain
  sourcemodel.pos   = sourcemodel.pos(sel,:);
  sourcemodel.dim   = [sum(xsel) sum(ysel) sum(zsel)];
end

fprintf('%d dipoles inside, %d dipoles outside brain\n', sum(sourcemodel.inside), sum(~sourcemodel.inside));

% apply the symmetry constraint, i.e. add a symmetric dipole for each location that was defined sofar
if ~isempty(cfg.symmetry)
  if size(sourcemodel.pos,2)>3
    % sanity check, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3119
    ft_warning('the construction of a symmetric dipole model requires to start with a Nx3 description of the dipole positions, discarding subsequent columns');
    sourcemodel.pos = sourcemodel.pos(:,1:3);
  end
  if strcmp(cfg.symmetry, 'x')
    reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    mirror = [1 1 1 -1 1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 -x1 y
  elseif strcmp(cfg.symmetry, 'y')
    reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    mirror = [1 1 1 1 -1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 -y
  elseif strcmp(cfg.symmetry, 'z')
    reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    mirror = [1 1 1 1 1 -1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 y1
  else
    ft_error('unrecognized symmetry constraint');
  end
  fprintf('each source describes two dipoles with symmetry along %s axis\n', cfg.symmetry);
  % expand the number of parameters from one (3) to two dipoles (6)
  sourcemodel.pos = sourcemodel.pos(:,expand) .* repmat(mirror, size(sourcemodel.pos,1), 1);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance sourcemodel
ft_postamble history    sourcemodel
ft_postamble savevar    sourcemodel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function for basedonmri method to determine the inside
% returns a boolean vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inside = getinside(pos, mask)

% it might be that the box with the points does not completely fit into the
% mask
dim = size(mask);
sel = find(pos(:,1)<1 |  pos(:,1)>dim(1) | ...
  pos(:,2)<1 |  pos(:,2)>dim(2) | ...
  pos(:,3)<1 |  pos(:,3)>dim(3));
if isempty(sel)
  % use the efficient implementation
  inside = mask(sub2ind(dim, pos(:,1), pos(:,2), pos(:,3)));
else
  % only loop over the points that can be dealt with
  inside = false(size(pos,1), 1);
  for i=setdiff(1:size(pos,1), sel(:)')
    inside(i) = mask(pos(i,1), pos(i,2), pos(i,3));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to return the fieldnames of the boolean fields in a
% segmentation, should work both for volumetric and for source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn = booleanfields(mri)

fn = fieldnames(mri);
isboolean = false(1,numel(fn));
for i=1:numel(fn)
  if islogical(mri.(fn{i})) && isequal(numel(mri.(fn{i})),prod(mri.dim))
    isboolean(i) = true;
  end
end
fn  = fn(isboolean);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute the centroids of the elements of a volumetric mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function centr = compute_centroids(headmodel)
centr = [];
% the FEM model should have tetraheders or hexaheders
if isfield(headmodel, 'tet')
  numtet = size(headmodel.tet, 1);
  centr.pos = zeros(numtet, 3);
  for i=1:numtet
    % compute the mean of the 4 corner points of the tetraheder
    centr.pos(i,:) = mean(headmodel.pos(hm.tet(i,:),:), 1);
  end
elseif isfield(headmodel, 'hex')
  numhex = size(headmodel.hex, 1);
  centr.pos = zeros(numhex, 3);
  for i=1:numhex
    % compute the mean of the 8 corner points of the hexaheder
    centr.pos(i,:) = mean(headmodel.pos(headmodel.hex(i,:),:), 1);
  end
else
  ft_error('the headmodel does not contain tetraheders or hexaheders');
end

% copy the specified fields, fields that are specified but not present will be silently ignored
centr = copyfields(headmodel, centr, {'tissue', 'tissuelabel', 'unit', 'coordsys'});
