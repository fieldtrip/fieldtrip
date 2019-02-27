function [sourcemodel, cfg] = ft_prepare_sourcemodel(cfg, headmodel, sens)

% FT_PREPARE_SOURCEMODEL constructs a source model, for example a 3D grid or a
% cortical sheet. The source model that can be used for source reconstruction,
% beamformer scanning, linear estimation and MEG interpolation.
%
% Use as
%   sourcemodel = ft_prepare_sourcemodel(cfg)
% where the configuration structure contains the details on how the source
% model should be constructed.
%
% A source model can be constructed based on
%   - regular 3D grid with explicit specification
%   - regular 3D grid with specification of the resolution
%   - regular 3D grid, based on segmented MRI, restricted to gray matter
%   - regular 3D grid, based on a warped template grid, based on the MNI brain
%   - surface mesh based on the brain surface from the volume conduction model
%   - surface mesh based on the head surface from an external file
%   - cortical sheet that was created in MNE or Freesurfer
%   - using user-supplied source positions, which can be regular or irregular
% The approach that will be used depends on the configuration options that
% you specify.
%
% Configuration options for generating a regular 3D grid
%   cfg.sourcemodel.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.sourcemodel.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.sourcemodel.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.sourcemodel.resolution = number (e.g. 1 cm) for automatic grid generation
%
% Configuration options for predefined source positions
%   cfg.sourcemodel.pos        = N*3 matrix with position of each source
%   cfg.sourcemodel.inside     = N*1 vector with boolean value whether position is inside brain (optional)
%   cfg.sourcemodel.dim        = [Nx Ny Nz] vector with dimensions in case of 3D grid (optional)
%
% The following fields are not used in this function, but will be copied along to the output
%   cfg.sourcemodel.leadfield
%   cfg.sourcemodel.filter or alternatively cfg.sourcemodel.avg.filter
%   cfg.sourcemodel.subspace
%   cfg.sourcemodel.lbex
%
% Configuration options for a warped MNI grid
%   cfg.mri                    = can be filename or MRI structure, containing the individual anatomy
%   cfg.sourcemodel.warpmni    = 'yes'
%   cfg.sourcemodel.resolution = number (e.g. 6) of the resolution of the
%                                template MNI grid, defined in mm
%   cfg.sourcemodel.template   = specification of a template grid (grid structure), or a
%                                filename of a template grid (defined in MNI space),
%                                either cfg.sourcemodel.resolution or cfg.sourcemodel.template needs
%                                to be defined. If both are defined cfg.sourcemodel.template prevails
%   cfg.sourcemodel.nonlinear  = 'no' (or 'yes'), use non-linear normalization
%
% Configuration options for cortex segmentation, i.e. for placing dipoles in grey matter
%   cfg.mri           = can be filename, MRI structure or segmented MRI structure
%   cfg.threshold     = 0.1, relative to the maximum value in the segmentation
%   cfg.smooth        = 5, smoothing in voxels
%
% Configuration options for reading a cortical sheet from file
%   cfg.headshape     = string, should be a *.fif file
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad          = structure with gradiometer definition or filename, see FT_READ_SENS
%
% The headmodel or volume conduction model can be specified as
%   cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%
% Other configuration options
%   cfg.sourcemodel.unit    = string, can be 'mm', 'cm', 'm' (default is automatic)
%   cfg.sourcemodel.tight   = 'yes' or 'no' (default is automatic)
%   cfg.inwardshift  = number, how much should the innermost surface be moved inward to constrain
%                      sources to be considered inside the source compartment (default = 0)
%   cfg.moveinward   = number, move dipoles inward to ensure a certain distance to the innermost
%                      surface of the source compartment (default = 0)
%   cfg.spherify     = 'yes' or 'no', scale the source model so that it fits inside a sperical
%                      volume conduction model (default = 'no')
%   cfg.symmetry     = 'x', 'y' or 'z' symmetry for two dipoles, can be empty (default = [])
%   cfg.headshape    = a filename for the headshape, a structure containing a single surface,
%                      or a Nx3 matrix with headshape surface points (default = [])
%   cfg.spmversion   = string, 'spm2', 'spm8', 'spm12' (default = 'spm8')
%
% See also FT_PREPARE_LEADFIELD, FT_PREPARE_HEADMODEL, FT_SOURCEANALYSIS,
% FT_DIPOLEFITTING, FT_MEGREALIGN

% Copyright (C) 2004-2013, Robert Oostenveld
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
cfg = ft_checkconfig(cfg, 'deprecated', 'mriunits');
cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'grid',    'sourcemodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed', {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed', {'optofile', 'opto'});

% put the low-level options pertaining to the sourcemodel in their own field
cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'}); % this is moved to cfg.sourcemodel.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.sourcemodel.unit by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'createsubcfg', {'sourcemodel'});

% set the defaults
cfg.moveinward        = ft_getopt(cfg, 'moveinward'); % the default is automatic and depends on a triangulation being present
cfg.spherify          = ft_getopt(cfg, 'spherify', 'no');
cfg.headshape         = ft_getopt(cfg, 'headshape');
cfg.symmetry          = ft_getopt(cfg, 'symmetry');
cfg.spmversion        = ft_getopt(cfg, 'spmversion', 'spm8');
cfg.sourcemodel       = ft_getopt(cfg, 'sourcemodel');
cfg.sourcemodel.unit  = ft_getopt(cfg.sourcemodel, 'unit', 'auto');

% this code expects the inside to be represented as a logical array
if isfield(cfg, 'sourcemodel')
  cfg.sourcemodel = ft_checkconfig(cfg.sourcemodel, 'renamed',  {'pnt' 'pos'});
  if isfield(cfg.sourcemodel, 'template')
    cfg.sourcemodel.template = ft_checkconfig(cfg.sourcemodel.template, 'renamed',  {'pnt' 'pos'});
  end
end
cfg = ft_checkconfig(cfg, 'inside2logical', 'yes');

if ~isfield(cfg, 'headmodel') && nargin>1
  % put it in the configuration structure
  % this is for backward compatibility, 13 Januari 2011
  cfg.headmodel = headmodel;
end

if ~isfield(cfg, 'grad') && ~isfield(cfg, 'elec') && nargin>2
  % put it in the configuration structure
  % this is for backward compatibility, 13 Januari 2011
  cfg.grad = sens;
end

if isfield(cfg.sourcemodel, 'resolution') && isfield(cfg.sourcemodel, 'xgrid') && ~ischar(cfg.sourcemodel.xgrid)
  ft_error('You cannot specify cfg.sourcemodel.resolution and an explicit cfg.sourcemodel.xgrid simultaneously');
end
if isfield(cfg.sourcemodel, 'resolution') && isfield(cfg.sourcemodel, 'ygrid') && ~ischar(cfg.sourcemodel.ygrid)
  ft_error('You cannot specify cfg.sourcemodel.resolution and an explicit cfg.sourcemodel.ygrid simultaneously');
end
if isfield(cfg.sourcemodel, 'resolution') && isfield(cfg.sourcemodel, 'zgrid') && ~ischar(cfg.sourcemodel.zgrid)
  ft_error('You cannot specify cfg.sourcemodel.resolution and an explicit cfg.sourcemodel.zgrid simultaneously');
end

% the source model can be constructed in a number of ways
basedongrid       = isfield(cfg.sourcemodel, 'xgrid') && ~ischar(cfg.sourcemodel.xgrid);                              % regular 3D grid with explicit specification
basedonpos        = isfield(cfg.sourcemodel, 'pos');                                                                  % using user-supplied positions, which can be regular or irregular
basedonshape      = ~isempty(cfg.headshape);                                                                          % surface mesh based on inward shifted head surface from external file
basedonmri        = isfield(cfg, 'mri') && ~(isfield(cfg.sourcemodel, 'warpmni') && istrue(cfg.sourcemodel.warpmni)); % regular 3D grid, based on segmented MRI, restricted to gray matter
basedonmni        = isfield(cfg, 'mri') &&  (isfield(cfg.sourcemodel, 'warpmni') && istrue(cfg.sourcemodel.warpmni)); % regular 3D grid, based on warped MNI template
basedonvol        = false;                                                                                            % surface mesh based on inward shifted brain surface from volume conductor
basedoncortex     = isfield(cfg, 'headshape') && (iscell(cfg.headshape) || any(ft_filetype(cfg.headshape, {'neuromag_fif', 'freesurfer_triangle_binary', 'caret_surf', 'gifti'}))); % cortical sheet from external software such as Caret or FreeSurfer, can also be two separate hemispheres
basedonresolution = isfield(cfg.sourcemodel, 'resolution') && ~basedonmri && ~basedonmni;                             % regular 3D grid with specification of the resolution

if basedonshape && basedoncortex
  % treating it as cortical sheet has preference
  basedonshape = false;
end

if basedongrid && basedonpos
  % fall back to default behaviour, in which the pos overrides the grid
  basedongrid = false;
end

if ~any([basedonresolution basedongrid basedonpos basedonshape basedonmri basedoncortex basedonmni]) && ~isempty(cfg.headmodel)
  % fall back to default behaviour, which is to create a surface grid (e.g. used in FT_MEGREALIGN)
  basedonvol = true;
end

% these are mutually exclusive, but printing all requested methods here
% facilitates debugging of weird configs. Also specify the defaults here to
% keep the overview
if basedonresolution
  fprintf('creating sourcemodel based on automatic 3D grid with specified resolution\n');
  cfg.sourcemodel.xgrid  = ft_getopt(cfg.sourcemodel, 'xgrid',  'auto');
  cfg.sourcemodel.ygrid  = ft_getopt(cfg.sourcemodel, 'ygrid',  'auto');
  cfg.sourcemodel.zgrid  = ft_getopt(cfg.sourcemodel, 'zgrid',  'auto');
  cfg.inwardshift = ft_getopt(cfg,      'inwardshift', 0); % in this case for inside detection, FIXME move to cfg.sourcemodel
  cfg.sourcemodel.tight  = ft_getopt(cfg.sourcemodel, 'tight',   'yes');
end

if basedongrid
  fprintf('creating sourcemodel based on user specified 3D grid\n');
  cfg.inwardshift = ft_getopt(cfg,      'inwardshift', 0); % in this case for inside detection, FIXME move to cfg.sourcemodel
  cfg.sourcemodel.tight  = ft_getopt(cfg.sourcemodel, 'tight',   'yes');
end

if basedonpos
  fprintf('creating sourcemodel based on user specified dipole positions\n');
  cfg.inwardshift = ft_getopt(cfg,      'inwardshift', 0); %in this case for inside detection, FIXME move to cfg.sourcemodel
  cfg.sourcemodel.tight  = ft_getopt(cfg.sourcemodel, 'tight',    'no');
end

if basedonshape
  fprintf('creating sourcemodel based on inward-shifted head shape\n');
  cfg.inwardshift = ft_getopt(cfg,      'inwardshift',  0); % in this case for inside detection, FIXME move to cfg.sourcemodel
  cfg.spheremesh  = ft_getopt(cfg,      'spheremesh', 642); % FIXME move spheremesh to cfg.sourcemodel
  cfg.sourcemodel.tight  = ft_getopt(cfg.sourcemodel, 'tight',    'yes');
end

if basedoncortex
  cfg.sourcemodel.tight  = ft_getopt(cfg.sourcemodel, 'tight', 'yes');
end

if basedonmri
  fprintf('creating sourcemodel based on an anatomical volume\n');
  cfg.threshold   = ft_getopt(cfg,      'threshold', 0.1);    % relative
  cfg.smooth      = ft_getopt(cfg,      'smooth',      5);    % in voxels
  cfg.sourcemodel.tight  = ft_getopt(cfg.sourcemodel, 'tight',   'yes');
end

if basedonvol
  fprintf('creating sourcemodel based on inward-shifted brain surface from volume conductor model\n');
  cfg.inwardshift = ft_getopt(cfg,      'inwardshift',   0); % in this case for inside detection, FIXME move to cfg.sourcemodel
  cfg.spheremesh  = ft_getopt(cfg,      'spheremesh',  642); % FIXME move spheremesh to cfg.sourcemodel
  cfg.sourcemodel.tight  = ft_getopt(cfg.sourcemodel, 'tight',      'no');
end

if basedonmni
  cfg.sourcemodel.tight     = ft_getopt(cfg.sourcemodel, 'tight',       'no');
  cfg.sourcemodel.nonlinear = ft_getopt(cfg.sourcemodel, 'nonlinear',   'no');
end

% these are mutually exclusive
if sum([basedonresolution basedongrid basedonpos basedonshape basedonmri basedonvol basedoncortex basedonmni])~=1
  ft_error('incorrect cfg specification for constructing a sourcemodel');
end

if (isfield(cfg, 'smooth') && ~strcmp(cfg.smooth, 'no')) || basedonmni
  % check that the preferred SPM version is on the path
  ft_hastoolbox(cfg.spmversion, 1);
end

% start with an empty structure
sourcemodel = [];

% get the volume conduction model
try
  headmodel = ft_fetch_headmodel(cfg);
catch
  headmodel = [];
end

% get the gradiometer or electrode definition
try
  sens = ft_fetch_sens(cfg);
catch
  sens = [];
end

if strcmp(cfg.sourcemodel.unit, 'auto')
  if isfield(cfg.sourcemodel, 'pos') && size(cfg.sourcemodel.pos,1)>10
    % estimate the units based on the existing source positions
    cfg.sourcemodel = rmfield(cfg.sourcemodel, 'unit'); % remove 'auto' and have ft_determine_units determine it properly
    cfg.sourcemodel = ft_determine_units(cfg.sourcemodel);
  elseif ~isempty(sens)
    % copy the units from the sensor array
    cfg.sourcemodel.unit = sens.unit;
  elseif ~isempty(headmodel)
    % copy the units from the volume conduction model
    cfg.sourcemodel.unit = headmodel.unit;
  else
    ft_warning('assuming "cm" as default source units');
    cfg.sourcemodel.unit = 'cm';
  end
end

% convert the sensor array to the desired units for the source model
if ~isempty(sens)
  sens = ft_convert_units(sens, cfg.sourcemodel.unit);
end

% convert the head model to the desired units for the source model
if ~isempty(headmodel)
  headmodel = ft_convert_units(headmodel, cfg.sourcemodel.unit);
end

if basedonresolution
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

  fprintf('creating 3D grid with %g %s resolution\n', cfg.sourcemodel.resolution, cfg.sourcemodel.unit);

  % round the bounding box limits to the nearest cm
  switch cfg.sourcemodel.unit
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

  if ischar(cfg.sourcemodel.xgrid) && strcmp(cfg.sourcemodel.xgrid, 'auto')
    sourcemodel.xgrid = minpos(1):cfg.sourcemodel.resolution:maxpos(1);
  end
  if ischar(cfg.sourcemodel.ygrid) && strcmp(cfg.sourcemodel.ygrid, 'auto')
    sourcemodel.ygrid = minpos(2):cfg.sourcemodel.resolution:maxpos(2);
  end
  if ischar(cfg.sourcemodel.zgrid) && strcmp(cfg.sourcemodel.zgrid, 'auto')
    sourcemodel.zgrid = minpos(3):cfg.sourcemodel.resolution:maxpos(3);
  end
  sourcemodel.dim   = [length(sourcemodel.xgrid) length(sourcemodel.ygrid) length(sourcemodel.zgrid)];
  [X, Y, Z]  = ndgrid(sourcemodel.xgrid, sourcemodel.ygrid, sourcemodel.zgrid);
  sourcemodel.pos   = [X(:) Y(:) Z(:)];
  sourcemodel.unit  = cfg.sourcemodel.unit;
  fprintf('initial 3D grid dimensions are [%d %d %d]\n', sourcemodel.dim(1), sourcemodel.dim(2), sourcemodel.dim(3));
end

if basedongrid
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % a detailed xgrid/ygrid/zgrid has been specified, the other details
  % still need to be determined
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sourcemodel.xgrid = cfg.sourcemodel.xgrid;
  sourcemodel.ygrid = cfg.sourcemodel.ygrid;
  sourcemodel.zgrid = cfg.sourcemodel.zgrid;
  sourcemodel.dim   = [length(sourcemodel.xgrid) length(sourcemodel.ygrid) length(sourcemodel.zgrid)];
  [X, Y, Z]  = ndgrid(sourcemodel.xgrid, sourcemodel.ygrid, sourcemodel.zgrid);
  sourcemodel.pos   = [X(:) Y(:) Z(:)];
  sourcemodel.unit  = cfg.sourcemodel.unit;
end

if basedonpos
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % source positions are already specified in the configuration, reuse as much of the
  % prespecified model as possible (but only known objects)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sourcemodel = keepfields(cfg.sourcemodel, {'pos', 'unit', 'xgrid', 'ygrid', 'zgrid', 'mom', 'tri', 'dim', 'transform', 'inside', 'lbex', 'subspace', 'leadfield', 'filter', 'label', 'leadfielddimord'});
end

if basedonmri
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

  if ~isfield(cfg.sourcemodel, 'resolution')
    switch cfg.sourcemodel.unit
      case 'mm'
        cfg.sourcemodel.resolution = 10;
      case 'cm'
        cfg.sourcemodel.resolution = 1;
      case 'dm'
        cfg.sourcemodel.resolution = 0.1;
      case 'm'
        cfg.sourcemodel.resolution = 0.01;
    end
  end

  issegmentation = false;
  if isfield(mri, 'gray')
    % this is not a boolean segmentation, but based on tissue probability
    % maps, being the original implementation here.
    dat = double(mri.gray);

    % apply a smoothing of a certain amount of voxels
    if ~strcmp(cfg.smooth, 'no');
      dat = volumesmooth(dat, cfg.smooth, 'MRI gray matter');
    end

  elseif isfield(mri, 'anatomy')
    % this could be a tpm stored on disk, i.e. the result of
    % ft_volumesegment. Reading it in always leads to the field 'anatomy'.
    % Note this could be any anatomical mask
    dat = double(mri.anatomy);

    % apply a smoothing of a certain amount of voxels
    if ~strcmp(cfg.smooth, 'no');
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
  scale = ft_scalingfactor(cfg.sourcemodel.unit, mri.unit);

  ind                 = find(head(:));
  fprintf('%d from %d voxels in the segmentation are marked as ''inside'' (%.0f%%)\n', length(ind), numel(head), 100*length(ind)/numel(head));
  [X,Y,Z]             = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));  % create the grid in MRI-coordinates
  posmri              = [X(ind) Y(ind) Z(ind)];                            % take only the inside voxels
  poshead             = ft_warp_apply(mri.transform, posmri);                 % transform to head coordinates
  resolution          = cfg.sourcemodel.resolution*scale;                                        % source and mri can be expressed in different units (e.g. cm and mm)
  xgrid               = floor(min(poshead(:,1))):resolution:ceil(max(poshead(:,1)));      % create the grid in head-coordinates
  ygrid               = floor(min(poshead(:,2))):resolution:ceil(max(poshead(:,2)));      % with 'consistent' x,y,z definitions
  zgrid               = floor(min(poshead(:,3))):resolution:ceil(max(poshead(:,3)));
  [X,Y,Z]             = ndgrid(xgrid,ygrid,zgrid);
  pos2head            = [X(:) Y(:) Z(:)];
  pos2mri             = ft_warp_apply(inv(mri.transform), pos2head);        % transform to MRI voxel coordinates
  pos2mri             = round(pos2mri);
  inside              = getinside(pos2mri, head);                           % use helper subfunction

  sourcemodel.pos            = pos2head/scale;                                     % convert to source units
  sourcemodel.xgrid          = xgrid/scale;                                        % convert to source units
  sourcemodel.ygrid          = ygrid/scale;                                        % convert to source units
  sourcemodel.zgrid          = zgrid/scale;                                        % convert to source units
  sourcemodel.dim            = [length(sourcemodel.xgrid) length(sourcemodel.ygrid) length(sourcemodel.zgrid)];
  sourcemodel.inside         = inside(:);
  sourcemodel.unit           = cfg.sourcemodel.unit;

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
end

if basedoncortex
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
  shape     = ft_convert_units(shape, cfg.sourcemodel.unit);
  % return both the vertices and triangles from the cortical sheet
  sourcemodel.pos  = shape.pos;
  sourcemodel.tri  = shape.tri;
  sourcemodel.unit = shape.unit;
end

if basedonshape
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the headshape  to make a superficial dipole layer (e.g.
  % for megrealign). Assume that all points are inside the volume.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % get the surface describing the head shape
  if isstruct(cfg.headshape) && isfield(cfg.headshape, 'hex')
    cfg.headshape = fixpos(cfg.headshape);
    fprintf('extracting surface from hexahedral mesh\n');
    headshape = mesh2edge(cfg.headshape);
    headshape = poly2tri(headshape);
  elseif isstruct(cfg.headshape) && isfield(cfg.headshape, 'tet')
    cfg.headshape = fixpos(cfg.headshape);
    fprintf('extracting surface from tetrahedral mesh\n');
    headshape = mesh2edge(cfg.headshape);
  elseif isstruct(cfg.headshape) && isfield(cfg.headshape, 'tri')
    cfg.headshape = fixpos(cfg.headshape);
    headshape = cfg.headshape;
  elseif isstruct(cfg.headshape) && isfield(cfg.headshape, 'pos')
    cfg.headshape = fixpos(cfg.headshape);
    headshape = cfg.headshape;
  elseif isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
    cfg.headshape = fixpos(cfg.headshape);
    headshape = cfg.headshape;
  elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
    % use the headshape points specified in the configuration
    cfg.headshape = fixpos(cfg.headshape);
    headshape = cfg.headshape;
  elseif ischar(cfg.headshape)
    % read the headshape from file
    headshape = ft_read_headshape(cfg.headshape);
  else
    ft_error('cfg.headshape is not specified correctly')
  end
  % ensure that the headshape is in the same units as the source model
  headshape = ft_convert_units(headshape, cfg.sourcemodel.unit);
  if ~isfield(headshape, 'tri') && ~isfield(headshape, 'poly')
    % generate a closed triangulation from the surface points
    headshape.pos = unique(headshape.pos, 'rows');
    headshape.tri = projecttri(headshape.pos);
  end
  % note that cfg.inwardshift should be expressed in the units consistent with cfg.sourcemodel.unit
  sourcemodel.pos     = headsurface([], [], 'headshape', headshape, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  sourcemodel.tri     = headshape.tri;
  sourcemodel.unit    = headshape.unit;
  sourcemodel.inside  = true(size(sourcemodel.pos,1),1);
end

if basedonvol
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the volume conduction model to make a superficial dipole layer (e.g.
  % for megrealign). Assume that all points are inside the volume.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % please note that cfg.inwardshift should be expressed in the units consistent with cfg.sourcemodel.unit
  sourcemodel.pos     = headsurface(headmodel, sens, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  sourcemodel.unit    = cfg.sourcemodel.unit;
  sourcemodel.inside  = true(size(sourcemodel.pos,1),1);
end

if basedonmni
  if ~isfield(cfg.sourcemodel, 'template') && ~isfield(cfg.sourcemodel, 'resolution')
    ft_error('you either need to specify the filename of a template grid in cfg.sourcemodel.template, or a resolution in cfg.sourcemodel.resolution');
  elseif isfield(cfg.sourcemodel, 'template')
    % let the template filename prevail
    fname = cfg.sourcemodel.template;
  elseif isfield(cfg.sourcemodel, 'resolution') && cfg.sourcemodel.resolution==round(cfg.sourcemodel.resolution)
    % use one of the templates that are in Fieldtrip, this requires a
    % resolution
    fname = ['standard_sourcemodel3d',num2str(cfg.sourcemodel.resolution),'mm.mat'];
  elseif isfield(cfg.sourcemodel, 'resolution') && cfg.sourcemodel.resolution~=round(cfg.sourcemodel.resolution)
    fname = ['standard_sourcemodel3d',num2str(floor(cfg.sourcemodel.resolution)),'point',num2str(10*(cfg.sourcemodel.resolution-floor(cfg.sourcemodel.resolution))),'mm.mat'];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % check whether the mni template grid exists for the specified resolution
  % if not create it: FIXME (this needs to be done still)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % get the mri
  if ischar(cfg.mri)
    if ~exist(fname, 'file')
      ft_error('the MNI template grid based on the specified resolution does not exist');
    end
    mri = ft_read_mri(cfg.mri);
  else
    mri = cfg.mri;
  end

  % get the template grid
  if ischar(fname)
    mnigrid = load(fname, 'sourcemodel');
    mnigrid = mnigrid.sourcemodel;
  else
    mnigrid = cfg.sourcemodel.template;
  end

  % ensure these to have units in mm, the conversion of the source model is done further down
  mri     = ft_convert_units(mri,     'mm');
  mnigrid = ft_convert_units(mnigrid, 'mm');

  % ensure that it is specified with logical inside
  mnigrid = fixinside(mnigrid);

  % spatial normalisation of mri and construction of subject specific sourcemodel positions
  tmpcfg           = keepfields(cfg, {'spmversion', 'spmmethod'});
  tmpcfg.nonlinear = cfg.sourcemodel.nonlinear;
  if isfield(cfg.sourcemodel, 'templatemri')
    tmpcfg.template = cfg.sourcemodel.templatemri;
  end
  normalise = ft_volumenormalise(tmpcfg, mri);

  if ~isfield(normalise, 'params') && ~isfield(normalise, 'initial')
    fprintf('applying an inverse warp based on a linear transformation only\n');
    sourcemodel.pos = ft_warp_apply(inv(normalise.cfg.final), mnigrid.pos);
  else
    sourcemodel.pos = ft_warp_apply(inv(normalise.initial), ft_warp_apply(normalise.params, mnigrid.pos, 'sn2individual'));
  end
  if isfield(mnigrid, 'dim')
    sourcemodel.dim   = mnigrid.dim;
  end
  if isfield(mnigrid, 'tri')
    sourcemodel.tri   = mnigrid.tri;
  end
  sourcemodel.unit    = mnigrid.unit;
  sourcemodel.inside  = mnigrid.inside;
  sourcemodel.params  = normalise.params;
  sourcemodel.initial = normalise.initial;
  if ft_datatype(mnigrid, 'parcellation')
    % copy the boolean fields over
    sourcemodel = copyfields(mnigrid, sourcemodel, booleanfields(mnigrid));
  end

end

% in most cases the source model will already be in the desired units, but e.g. for "basedonmni" it will be in 'mm'
% convert to the requested units
sourcemodel = ft_convert_units(sourcemodel, cfg.sourcemodel.unit);

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

% determine the dipole locations that are inside the source compartment of the
% volume conduction model, i.e. inside the brain
if ~isfield(sourcemodel, 'inside')
  sourcemodel.inside = ft_inside_headmodel(sourcemodel.pos, headmodel, 'grad', sens, 'headshape', cfg.headshape, 'inwardshift', cfg.inwardshift); % this returns a boolean vector
end

if strcmp(cfg.sourcemodel.tight, 'yes')
  fprintf('%d dipoles inside, %d dipoles outside brain\n', sum(sourcemodel.inside), sum(~sourcemodel.inside));
  fprintf('making tight grid\n');
  xmin = min(sourcemodel.pos(sourcemodel.inside,1));
  ymin = min(sourcemodel.pos(sourcemodel.inside,2));
  zmin = min(sourcemodel.pos(sourcemodel.inside,3));
  xmax = max(sourcemodel.pos(sourcemodel.inside,1));
  ymax = max(sourcemodel.pos(sourcemodel.inside,2));
  zmax = max(sourcemodel.pos(sourcemodel.inside,3));
  xmin_indx = find(sourcemodel.xgrid==xmin);
  ymin_indx = find(sourcemodel.ygrid==ymin);
  zmin_indx = find(sourcemodel.zgrid==zmin);
  xmax_indx = find(sourcemodel.xgrid==xmax);
  ymax_indx = find(sourcemodel.ygrid==ymax);
  zmax_indx = find(sourcemodel.zgrid==zmax);
  sel =       (sourcemodel.pos(:,1)>=xmin & sourcemodel.pos(:,1)<=xmax); % select all grid positions inside the tight box
  sel = sel & (sourcemodel.pos(:,2)>=ymin & sourcemodel.pos(:,2)<=ymax); % select all grid positions inside the tight box
  sel = sel & (sourcemodel.pos(:,3)>=zmin & sourcemodel.pos(:,3)<=zmax); % select all grid positions inside the tight box
  % update the grid locations that are marked as inside the brain
  sourcemodel.pos   = sourcemodel.pos(sel,:);
  % update the boolean fields, this requires the original dim
  fn = booleanfields(sourcemodel);
  for i=1:numel(fn)
    sourcemodel.(fn{i}) = sourcemodel.(fn{i})(sel);
  end
  sourcemodel.xgrid   = sourcemodel.xgrid(xmin_indx:xmax_indx);
  sourcemodel.ygrid   = sourcemodel.ygrid(ymin_indx:ymax_indx);
  sourcemodel.zgrid   = sourcemodel.zgrid(zmin_indx:zmax_indx);
  sourcemodel.dim     = [length(sourcemodel.xgrid) length(sourcemodel.ygrid) length(sourcemodel.zgrid)];
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

%--------------------------------------------------------------
% helper function for basedonmri method to determine the inside
% returns a boolean vector
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

%--------------------------------------------------------------------------
% helper function to return the fieldnames of the boolean fields in a
% segmentation, should work both for volumetric and for source
function fn = booleanfields(mri)

fn = fieldnames(mri);
isboolean = false(1,numel(fn));
for i=1:numel(fn)
  if islogical(mri.(fn{i})) && isequal(numel(mri.(fn{i})),prod(mri.dim))
    isboolean(i) = true;
  end
end
fn  = fn(isboolean);
