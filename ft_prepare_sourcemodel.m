function [grid, cfg] = ft_prepare_sourcemodel(cfg, vol, sens)

% FT_PREPARE_SOURCEMODEL helps to make a source model that can be
% used for source reconstruction, beamformer scanning, linear
% estimation and MEG interpolation.
%
% Use as
%   [grid, cfg] = prepare_dipole_grid(cfg)
% where the configuration structure contains the details on how the source
% model should be constructed.
%
% A source model can be constructed based on
%   - regular 3D grid with explicit specification
%   - regular 3D grid with specification of the resolution
%   - regular 3D grid, based on segmented MRI, restricted to gray matter
%   - surface grid based on brain surface from volume conductor
%   - surface grid based on head surface from external file
%   - using user-supplied grid positions, which can be regular or irregular
%   - cortical sheet that was created in MNE or Freesurfer
% The approach that will be used depends on the configuration options that
% you specify.
%
% Configuration options for generating a regular 3-D grid
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
%
% Configuration options for a predefined grid
%   cfg.grid.pos        = N*3 matrix with position of each source
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%   cfg.grid.inside     = vector with indices of the sources inside the brain (optional)
%   cfg.grid.outside    = vector with indices of the sources outside the brain (optional)
% The following fields are not used in this function, but will be copied along to the output
%   cfg.grid.leadfield
%   cfg.grid.filter or alternatively cfg.grid.avg.filter
%   cfg.grid.subspace
%   cfg.grid.lbex
%
% Configuration options for cortex segmentation, i.e. for placing dipoles in grey matter
%   cfg.mri           = can be filename, MRI structure or segmented MRI structure
%   cfg.sourceunits   = 'mm' or 'cm' (default is 'cm')
%   cfg.mriunits      = 'mm' or 'cm' (default is 'mm')
%   cfg.threshold     = 0.1, relative to the maximum value in the segmentation
%   cfg.smooth        = 5, smoothing in voxels
%
% Configuration options for reading a cortical sheet from file
%   cfg.headshape     = string, should be a *.fif file
%
% Other configuration options
%   cfg.vol          = volume conduction model
%   cfg.grad         = gradiometer definition
%   cfg.elec         = electrode definition
%   cfg.grid.tight   = 'yes' or 'no' (default is automatic)
%   cfg.inwardshift  = depth of the bounding layer for the source space, relative to the head model surface (default = 0)
%   cfg.symmetry     = 'x', 'y' or 'z' symmetry for two dipoles, can be empty (default = [])
%   cfg.headshape    = a filename containing headshape, a structure containing a
%                      single triangulated boundary, or a Nx3 matrix with surface
%                      points
%
% See also FT_SOURCEANALYSIS, FT_MEGREALIGN, FT_DIPOLEFITTING, FT_PREPARE_LEADFIELD

% Searching through the code, it seems that the following cfg fields are being used
% cfg.grid
% cfg.mri
% cfg.headshape
% cfg.tightgrid
% cfg.symmetry
% cfg.smooth
% cfg.threshold
% cfg.spheremesh
% cfg.inwardshift
% cfg.mriunits
% cfg.sourceunits

% Copyright (C) 2004-2011, Robert Oostenveld
%
% %Log$

% set the defaults
if ~isfield(cfg, 'symmetry'),         cfg.symmetry    = [];       end
if ~isfield(cfg, 'grid'),             cfg.grid        = [];       end

if ~isfield(cfg, 'vol') && nargin>1
  % put it in the configuration structure
  % this is for backward compatibility, 13 Januari 2011
  cfg.vol = vol;
end

if ~isfield(cfg, 'sens') && nargin>2
  % put it in the configuration structure
  % this is for backward compatibility, 13 Januari 2011
  cfg.grad = sens;
end

% for backward compatibility
if isfield(cfg, 'tightgrid')
  cfg.grid.tight = cfg.tightgrid;
  cfg = rmfield(cfg, 'tightgrid');
end

if isfield(cfg.grid, 'resolution') && isfield(cfg.grid, 'xgrid') && ~ischar(cfg.grid.xgrid)
  error('You cannot specify cfg.grid.resolution and an explicit cfg.grid.xgrid simultaneously');
end
if isfield(cfg.grid, 'resolution') && isfield(cfg.grid, 'ygrid') && ~ischar(cfg.grid.ygrid)
  error('You cannot specify cfg.grid.resolution and an explicit cfg.grid.ygrid simultaneously');
end
if isfield(cfg.grid, 'resolution') && isfield(cfg.grid, 'zgrid') && ~ischar(cfg.grid.zgrid)
  error('You cannot specify cfg.grid.resolution and an explicit cfg.grid.zgrid simultaneously');
end

% a grid can be constructed based on a number of ways
basedonauto   = isfield(cfg.grid, 'resolution');                        % regular 3D grid with specification of the resolution
basedongrid   = isfield(cfg.grid, 'xgrid') && ~ischar(cfg.grid.xgrid);  % regular 3D grid with explicit specification
basedonpos    = isfield(cfg.grid, 'pos');                               % using user-supplied grid positions, which can be regular or irregular
basedonshape  = isfield(cfg, 'headshape') && ~isempty(cfg.headshape);   % surface grid based on inward shifted head surface from external file
basedonmri    = isfield(cfg, 'mri');                                    % regular 3D grid, based on segmented MRI, restricted to gray matter
basedonvol    = false;                                                  % surface grid based on inward shifted brain surface from volume conductor
basedoncortex = isfield(cfg, 'headshape') && (iscell(cfg.headshape) || ft_filetype(cfg.headshape, 'neuromag_fif') || ft_filetype(cfg.headshape, 'freesurfer_triangle_binary')); % cortical sheet from MNE or Freesurfer, also in case of multiple files/hemispheres

if basedonshape && basedoncortex
  % treating it as cortical sheet has preference
  basedonshape = false;
end

if basedongrid && basedonpos
  % fall back to default behaviour, in which the pos overrides the grid
  basedongrid = false;
end

if ~any([basedonauto basedongrid basedonpos basedonshape basedonmri basedoncortex]) && nargin>1 && ~isempty(cfg.vol)
  % fall back to default behaviour, which is to create a surface grid (e.g. used in MEGRELIGN)
  basedonvol = 1;
end

% these are mutually exclusive, but printing all requested methods here
% facilitates debugging of weird configs. Also specify the defaults here to
% keep the overview
if basedonauto
  fprintf('creating dipole grid based on automatic 3D grid with specified resolution\n');
  if ~isfield(cfg.grid, 'xgrid'),       cfg.grid.xgrid = 'auto';    end
  if ~isfield(cfg.grid, 'ygrid'),       cfg.grid.ygrid = 'auto';    end
  if ~isfield(cfg.grid, 'zgrid'),       cfg.grid.zgrid = 'auto';    end
  if ~isfield(cfg, 'inwardshift'),      cfg.inwardshift = 0;        end % in this case for inside detection, % FIXME move to cfg.grid
  if ~isfield(cfg.grid, 'tight'),       cfg.grid.tight = 'yes';     end
end

if basedongrid
  fprintf('creating dipole grid based on user specified 3D grid\n');
  if ~isfield(cfg, 'inwardshift'),      cfg.inwardshift = 0;        end % in this case for inside detection, % FIXME move to cfg.grid
  if ~isfield(cfg.grid, 'tight'),       cfg.grid.tight = 'no';      end
end

if basedonpos
  fprintf('creating dipole grid based on user specified dipole positions\n');
  if ~isfield(cfg, 'inwardshift'),      cfg.inwardshift = 0;        end % in this case for inside detection, % FIXME move to cfg.grid
  if ~isfield(cfg.grid, 'tight'),       cfg.grid.tight = 'no';      end
end

if basedonshape
  fprintf('creating dipole grid based on inward-shifted head shape\n');
  if ~isfield(cfg, 'spheremesh'),       cfg.spheremesh = 642;  end  % FIXME move spheremesh to cfg.grid
  if ~isfield(cfg.grid, 'tight'),       cfg.grid.tight = 'no'; end
end

if basedoncortex
  if ~isfield(cfg, 'sourceunits');      cfg.sourceunits = 'cm';     end
  if ~isfield(cfg.grid, 'tight'),       cfg.grid.tight = 'no'; end
end

if basedonmri
  fprintf('creating dipole grid based on grey matter from segmented MRI\n');
  if ~isfield(cfg, 'threshold'),        cfg.threshold = 0.1;        end % relative
  if ~isfield(cfg, 'smooth');           cfg.smooth    = 5;          end % in voxels
  if ~isfield(cfg, 'sourceunits');      cfg.sourceunits = 'cm';     end
  if ~isfield(cfg, 'mriunits');         cfg.mriunits = 'mm';        end
  if ~isfield(cfg.grid, 'tight'),       cfg.grid.tight = 'yes';     end
end

if basedonvol
  fprintf('creating dipole grid based on inward-shifted brain surface from volume conductor model\n');
  if ~isfield(cfg, 'spheremesh'),       cfg.spheremesh = 642;  end   % FIXME move to cfg.grid
  if ~isfield(cfg.grid, 'tight'),       cfg.grid.tight = 'no'; end
end

% these are mutually exclusive
if sum([basedonauto basedongrid basedonpos basedonshape basedonmri basedonvol basedoncortex])~=1
  error('incorrect cfg specification for constructing a dipole grid');
end

needspm = isfield(cfg, 'smooth') && ~strcmp(cfg.smooth, 'no');
if needspm
  % check if SPM is in path, do not add any of them to the path yet
  hasspm2 = ft_hastoolbox('SPM2');
  hasspm8 = ft_hastoolbox('SPM8');
  
  if ~hasspm2 && ~hasspm8
    % is neither is on the path, try adding SPM8
    try, hasspm8 = ft_hastoolbox('SPM8', 1); end
  end
  
  if ~hasspm8
    % if SPM8 could not be added, try adding SPM2
    try, ft_hastoolbox('SPM2', 1); end
  end
end

% start with an empty grid
grid = [];

% copy the volume conductor and sensor array out of the cfg
if isfield(cfg, 'vol')
  vol = cfg.vol;
else
  vol = [];
end

% these are mutually exclusive
if isfield(cfg, 'grad')
  sens = cfg.grad;
elseif isfield(cfg, 'elec')
  sens = cfg.elec;
else
  sens = [];
end

if basedonauto
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % construct a regular 3D grid that spans a box encompassing all electrode
  % or gradiometer coils, this will typically also cover the complete brain
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ischar(cfg.grid.xgrid) && strcmp(cfg.grid.xgrid, 'auto')
    grid.xgrid = floor(min(sens.pnt(:,1))):cfg.grid.resolution:ceil(max(sens.pnt(:,1)));
  end
  if ischar(cfg.grid.ygrid) && strcmp(cfg.grid.ygrid, 'auto')
    grid.ygrid = floor(min(sens.pnt(:,2))):cfg.grid.resolution:ceil(max(sens.pnt(:,2)));
  end
  if ischar(cfg.grid.zgrid) && strcmp(cfg.grid.zgrid, 'auto')
    grid.zgrid = floor(min(sens.pnt(:,3))):cfg.grid.resolution:ceil(max(sens.pnt(:,3)));
  end
  grid.dim   = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
  [X, Y, Z]  = ndgrid(grid.xgrid, grid.ygrid, grid.zgrid);
  grid.pos   = [X(:) Y(:) Z(:)];
end

if basedongrid
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % a detailled xgrid/ygrid/zgrid has been specified, the other details
  % still need to be determined
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  grid.xgrid = cfg.grid.xgrid;
  grid.ygrid = cfg.grid.ygrid;
  grid.zgrid = cfg.grid.zgrid;
  grid.dim   = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
  [X, Y, Z]  = ndgrid(grid.xgrid, grid.ygrid, grid.zgrid);
  grid.pos   = [X(:) Y(:) Z(:)];
end

if basedonpos
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % a grid is already specified in the configuration, reuse as much of the
  % prespecified grid as possible (but only known objects)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if isfield(cfg.grid, 'pos')
    grid.pos = cfg.grid.pos;
  end
  if isfield(cfg.grid, 'mom')
    grid.mom = cfg.grid.mom;
  end
  if isfield(cfg.grid, 'dim')
    grid.dim = cfg.grid.dim;
  end
  if isfield(cfg.grid, 'xgrid')
    % FIXME is it desirable to have this in the grid?
    grid.xgrid = cfg.grid.xgrid;
  end
  if isfield(cfg.grid, 'ygrid')
    % FIXME is it desirable to have this in the grid?
    grid.ygrid = cfg.grid.ygrid;
  end
  if isfield(cfg.grid, 'zgrid')
    % FIXME is it desirable to have this in the grid?
    grid.zgrid = cfg.grid.zgrid;
  end
  if isfield(cfg.grid, 'dim')
    grid.dim = cfg.grid.dim;
  elseif isfield(grid, 'xgrid') && isfield(grid, 'ygrid') && isfield(grid, 'zgrid')
    grid.dim = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
  end
  if isfield(cfg.grid, 'inside')
    grid.inside = cfg.grid.inside;
  end
  if isfield(cfg.grid, 'outside')
    grid.outside = cfg.grid.outside;
  end
  if isfield(cfg.grid, 'lbex')
    grid.lbex = cfg.grid.lbex;
  end
  if isfield(cfg.grid, 'subspace')
    grid.subspace = cfg.grid.subspace;
  end
  if isfield(cfg.grid, 'leadfield')
    grid.leadfield = cfg.grid.leadfield;
  end
  if isfield(cfg.grid, 'filter')
    grid.filter = cfg.grid.filter;
  end
  if isfield(cfg.grid, 'tight')
    grid.tight = cfg.grid.tight;
  end
  % this is not supported any more
  if isfield(cfg.grid, 'avg') && isfield(cfg.grid.avg, 'filter')
    error('please put your filters in cfg.grid instead of cfg.grid.avg');
  end
end

if basedonmri
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % construct a grid based on the segmented MRI that is provided in the
  % configuration, only voxels in gray matter will be used
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % convert the source/functional data into the same units as the anatomical MRI
  scale = 1;
  switch cfg.sourceunits
    case 'mm'
      scale = scale / 1000;
    case 'cm'
      scale = scale / 100;
    case 'dm'
      scale = scale / 10;
    case 'm'
      scale = scale / 1;
    otherwise
      error('unknown physical dimension in cfg.sourceunits');
  end
  switch cfg.mriunits
    case 'mm'
      scale = scale * 1000;
    case 'cm'
      scale = scale * 100;
    case 'dm'
      scale = scale * 10;
    case 'm'
      scale = scale * 1;
    otherwise
      error('unknown physical dimension in cfg.mriunits');
  end
  
  if ~isfield(cfg.grid, 'resolution')
    switch cfg.sourceunits
      case 'mm'
        cfg.grid.resolution = 10;
      case 'cm'
        cfg.grid.resolution = 1;
      case 'dm'
        cfg.grid.resolution = 0.1;
      case 'm'
        cfg.grid.resolution = 0.01;
    end
  end
  
  if ischar(cfg.mri)
    % read the segmentation from file
    mri      = ft_read_mri(cfg.mri);
    mri.gray = double(mri.anatomy);
  elseif isstruct(cfg.mri) && ~isfield(cfg.mri, 'gray')
    % looks like a segmentation that has already been loaded in memory
    mri      = cfg.mri;
    mri.gray = double(mri.anatomy);
  elseif isstruct(cfg.mri) && isfield(cfg.mri, 'gray')
    % looks like a complete segmentation from VOLUMESEGMENT
    mri      = cfg.mri;
    mri.gray = double(mri.gray);
  else
    error('cannot determine the format of the segmentation in cfg.mri');
  end
  
  % apply a smoothing of a certain amount of voxels
  if ~strcmp(cfg.smooth, 'no');
    fprintf('smoothing gray matter segmentation with %d voxels\n', cfg.smooth);
    spm_smooth(mri.gray, mri.gray, cfg.smooth);
  end
  
  % determine for each voxel whether it belongs to the cortex
  if isfield(cfg, 'threshold')
    fprintf('thresholding gray matter segmentation at a relative value of %f\n', cfg.threshold);
    head = mri.gray./max(mri.gray(:)) > cfg.threshold;
  else
    error('you must specify cfg.threshold for cortex segmentation');
  end
  
  ind                 = find(head(:));
  fprintf('%d from %d voxels in the segmentation are marked as cortex (%.0f%%)\n', length(ind), prod(size(head)), 100*length(ind)/prod(size(head)));
  [X,Y,Z]             = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));   % create the grid in MRI-coordinates
  posmri              = [X(ind) Y(ind) Z(ind) ones(length(ind),1)];         % take only the inside voxels
  poshead             = mri.transform * posmri';                            % transform to head coordinates
  poshead             = poshead(1:3,:)';
  posmri              = posmri(:,1:3);
  resolution          = cfg.grid.resolution*scale;                                        % source and mri can be expressed in different units (e.g. cm and mm)
  xgrid               = floor(min(poshead(:,1))):resolution:ceil(max(poshead(:,1)));      % create the grid in head-coordinates
  ygrid               = floor(min(poshead(:,2))):resolution:ceil(max(poshead(:,2)));      % with 'consistent' x,y,z definitions
  zgrid               = floor(min(poshead(:,3))):resolution:ceil(max(poshead(:,3)));
  [X,Y,Z]             = ndgrid(xgrid,ygrid,zgrid);
  pos2head            = [X(:) Y(:) Z(:) ones(length(X(:)),1)]';
  pos2mri             = mri.transform \ pos2head;                                         % transform to MRI-coordinates
  pos2mri             = round(pos2mri(1:3,:))';
  pos2head            = pos2head(1:3,:)';
  pos2mri             = pos2mri(:,1:3);
  % it might be that the box with the points does not completely fit into the MRI
  sel = find(pos2mri(:,1)<1 |  pos2mri(:,1)>size(head,1) | ...
    pos2mri(:,2)<1 |  pos2mri(:,2)>size(head,2) | ...
    pos2mri(:,3)<1 |  pos2mri(:,3)>size(head,3));
  if isempty(sel)
    % use the efficient implementation
    inside = head(sub2ind(mri.dim, pos2mri(:,1), pos2mri(:,2), pos2mri(:,3)));
  else
    % only loop over the points that can be dealt with
    inside = zeros(length(xgrid)*length(ygrid)*length(zgrid), 1);
    for i=setdiff(1:size(pos2mri,1), sel(:)')
      inside(i) = head(pos2mri(i,1), pos2mri(i,2), pos2mri(i,3));
    end
  end
  inside = find(inside);
  
  grid.pos            = pos2head/scale;                                     % convert to source units
  grid.xgrid          = xgrid/scale;                                        % convert to source units
  grid.ygrid          = ygrid/scale;                                        % convert to source units
  grid.zgrid          = zgrid/scale;                                        % convert to source units
  grid.dim            = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
  grid.inside         = inside(:);
  grid.outside        = setdiff(1:size(grid.pos,1),grid.inside)';
  
  fprintf('the regular 3D grid encompassing the cortex contains %d grid points\n', size(grid.pos,1));
  fprintf('%d grid points inside gray matter\n', length(grid.inside));
  fprintf('%d grid points outside gray matter\n', length(grid.outside));
end

if basedoncortex
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read it from a *.fif file that was created using Freesurfer and MNE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if iscell(cfg.headshape)
    % FIXME loop over all files, this should be two hemispheres
    keyboard
  end
  switch ft_filetype(cfg.headshape)
    case 'freesurfer_triangle_binary'
      % it contains a cortical sheet which was created by the Freesurfer software
      shape = ft_read_headshape(cfg.headshape, 'format', 'freesurfer_triangle_binary');
    case 'neuromag_fif'
      % fif files can contain a variety of objects
      % here we assume that it contains a cortical sheet which was created by the MNE software
      shape = ft_read_headshape(cfg.headshape, 'format', 'mne_source');
    otherwise
      % use autodetection
      shape = ft_read_headshape(cfg.headshape);
  end
  grid.pos = shape.pnt;
  grid.tri = shape.tri;
  grid = ft_convert_units(grid, cfg.sourceunits);
end

if basedonshape
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the headshape  to make a superficial dipole layer (e.g.
  % for megrealign). Assume that all points are inside the volume.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % get the surface describing the head shape
  if isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
    % use the headshape surface specified in the configuration
    headshape = cfg.headshape;
  elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
    % use the headshape points specified in the configuration
    headshape.pnt = cfg.headshape;
  elseif ischar(cfg.headshape)
    % read the headshape from file
    headshape = ft_read_headshape(cfg.headshape);
  else
    error('cfg.headshape is not specified correctly')
  end
  if ~isfield(headshape, 'tri')
    % generate a closed triangulation from the surface points
    headshape.pnt = unique(headshape.pnt, 'rows');
    headshape.tri = projecttri(headshape.pnt);
  end
  grid.pos = headsurface([], [], 'headshape', headshape, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  grid.inside  = 1:size(grid.pos,1);
  grid.outside = [];
end

if basedonvol
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the volume conduction model to make a superficial dipole layer (e.g.
  % for megrealign). Assume that all points are inside the volume.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  grid.pos = headsurface(vol, sens, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  grid.inside  = 1:size(grid.pos,1);
  grid.outside = [];
end

% FIXME use inside_vol instead of this replication of code
% determine the dipole locations inside the brain volume
if ~isfield(grid, 'inside') && ~isfield(grid, 'outside')
  if ft_voltype(vol, 'infinite') || ft_voltype(vol, 'halfspace')
    % an empty vol in combination with gradiometers indicates a magnetic dipole
    % in an infinite vacuum, i.e. all dipoles can be considered to be inside
    grid.inside = 1:size(grid.pos,1);
    grid.outside = [];
    outside = zeros(1,size(grid.pos,1));
    for i =1:size(grid.pos,1);
      invacuum = false;
      dip1 = grid.pos(i,:);
      % condition of dipoles falling in the non conductive halfspace
      invacuum = acos(dot(vol.ori,(dip1-vol.pnt)./norm(dip1-vol.pnt))) < pi/2;
      if invacuum
        outside(i) = 1;
      end
    end
    grid.outside = find(outside);
    grid.inside  = find(~outside);
  else
    if isfield(sens, 'ori') && isfield(sens, 'pnt') && isfield(sens, 'tra')
      % in case of MEG, make a triangulation of the outermost surface
      if isfield(cfg, 'headshape')
        % use the specified headshape to construct the bounding triangulation
        [pnt, tri] = headsurface(vol, sens, 'headshape', cfg.headshape, 'inwardshift', cfg.inwardshift, 'surface', 'skin');
      else
        % use the volume conductor model to construct the bounding triangulation
        [pnt, tri] = headsurface(vol, sens, 'inwardshift', cfg.inwardshift, 'surface', 'skin');
      end
    else
      % in case of EEG, make a triangulation of the innermost surface
      [pnt, tri] = headsurface(vol, sens, 'inwardshift', cfg.inwardshift, 'surface', 'brain');
    end
    % determine the dipole positions that are inside the triangulated surface
    tmp = bounding_mesh(grid.pos, pnt, tri);
    grid.inside  = find(tmp==1);
    grid.outside = find(tmp==0);
  end
elseif ~isfield(grid, 'inside')
  grid.inside = setdiff(1:size(grid.pos,1), grid.outside);
elseif ~isfield(grid, 'outside')
  grid.outside = setdiff(1:size(grid.pos,1), grid.inside);
end

if strcmp(cfg.grid.tight, 'yes')
  fprintf('%d dipoles inside, %d dipoles outside brain\n', length(grid.inside), length(grid.outside));
  fprintf('making tight grid\n');
  xmin = min(grid.pos(grid.inside,1));
  ymin = min(grid.pos(grid.inside,2));
  zmin = min(grid.pos(grid.inside,3));
  xmax = max(grid.pos(grid.inside,1));
  ymax = max(grid.pos(grid.inside,2));
  zmax = max(grid.pos(grid.inside,3));
  xmin_indx = find(grid.xgrid==xmin);
  ymin_indx = find(grid.ygrid==ymin);
  zmin_indx = find(grid.zgrid==zmin);
  xmax_indx = find(grid.xgrid==xmax);
  ymax_indx = find(grid.ygrid==ymax);
  zmax_indx = find(grid.zgrid==zmax);
  sel =       (grid.pos(:,1)>=xmin & grid.pos(:,1)<=xmax); % select all grid positions inside the tight box
  sel = sel & (grid.pos(:,2)>=ymin & grid.pos(:,2)<=ymax); % select all grid positions inside the tight box
  sel = sel & (grid.pos(:,3)>=zmin & grid.pos(:,3)<=zmax); % select all grid positions inside the tight box
  grid.pos   = grid.pos(sel,:);
  % update the inside and outside vector
  tmp = zeros(1,prod(grid.dim));
  tmp(grid.inside)  = 1;        % these are originally inside the brain
  tmp(grid.outside) = 0;        % these are originally outside the brain
  tmp               = tmp(sel); % within the tight box, these are inside the brain
  grid.inside  = find(tmp);
  grid.outside = find(~tmp);
  grid.xgrid   = grid.xgrid(xmin_indx:xmax_indx);
  grid.ygrid   = grid.ygrid(ymin_indx:ymax_indx);
  grid.zgrid   = grid.zgrid(zmin_indx:zmax_indx);
  grid.dim     = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
end
fprintf('%d dipoles inside, %d dipoles outside brain\n', length(grid.inside), length(grid.outside));

% apply the symmetry constraint, i.e. add a symmetric dipole for each location defined sofar
% set up the symmetry constraints
if ~isempty(cfg.symmetry)
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
    error('unrecognized symmetry constraint');
  end
  fprintf('each source describes two dipoles with symmetry along %s axis\n', cfg.symmetry);
  % expand the number of parameters from one (3) to two dipoles (6)
  grid.pos = grid.pos(:,expand) .* repmat(mirror, size(grid.pos,1), 1);
end

% FIXME should the cfg be added to the output grid?
