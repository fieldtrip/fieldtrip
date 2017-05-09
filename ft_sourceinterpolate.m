function [interp] = ft_sourceinterpolate(cfg, functional, anatomical)

% FT_SOURCEINTERPOLATE interpolates source activity or statistical maps onto the
% voxels or vertices of an anatomical description of the brain.  Both the functional
% and the anatomical data can either describe a volumetric 3D regular grid, a
% triangulated description of the cortical sheet or an arbitrary cloud of points.
%
% The functional data in the output data will be interpolated at the locations at
% which the anatomical data are defined. For example, if the anatomical data was
% volumetric, the output data is a volume-structure, containing the resliced source
% and the anatomical volume that can be visualized using FT_SOURCEPLOT or written to
% file using FT_SOURCEWRITE.
%
% The following scenarios are possible:
%
% - Both functional data and anatomical data are defined on 3D regular grids, for
%   example with a low-res grid for the functional data and a high-res grid for the
%   anatomy.
%
% - The functional data is defined on a 3D regular grid of source positions
%   and the anatomical data is defined on an irregular point cloud, which can be a
%   2D triangulated mesh.
%
% - The functional data is defined on an irregular point cloud, which can be a 2D
%   triangulated mesh, and the anatomical data is defined on a 3D regular grid.
%
% - Both the functional and the anatomical data are defined on an irregular
%   point cloud, which can be a 2D triangulated mesh.
%
% - The functional data is defined on a low resolution 2D triangulated mesh and the
%   anatomical data is defined on a high resolution mesh, where the low-res vertices
%   form a subset of the high-res vertices. This allows for mesh based interpolation.
%   The algorithm currently implemented is so-called 'smudging' as it is also applied
%   by the MNE-suite software.
%
% Use as
%   [interp] = ft_sourceinterpolate(cfg, source, anatomy)
%   [interp] = ft_sourceinterpolate(cfg, stat,   anatomy)
% where
%   source  is the output of FT_SOURCEANALYSIS
%   stat    is the output of FT_SOURCESTATISTICS
%   anatomy is the output of FT_READ_MRI or one of the FT_VOLUMExxx functions,
%           a cortical sheet that was read with FT_READ_HEADSHAPE, or a regular
%           3D grid created with FT_PREPARE_SOURCEMODEL.
% and cfg is a structure with any of the following fields
%   cfg.parameter     = string (or cell-array) of the parameter(s) to be interpolated
%   cfg.downsample    = integer number (default = 1, i.e. no downsampling)
%   cfg.interpmethod  = string, can be 'nearest', 'linear', 'cubic',  'spline', 'sphere_avg' or 'smudge' (default = 'linear for interpolating two 3D volumes, 'nearest' for all other cases)
%
% The supported interpolation methods are 'nearest', 'linear', 'cubic' or 'spline'
% for interpolating two 3D volumes onto each other. For all other cases the supported
% interpolation methods are 'nearest', 'sphere_avg' or 'smudge'.
%
% The functional and anatomical data should be expressed in the same
% coordinate sytem, i.e. either both in MEG headcoordinates (NAS/LPA/RPA)
% or both in SPM coordinates (AC/PC).
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_READ_MRI, FT_SOURCEANALYSIS, FT_SOURCESTATISTICS,
% FT_READ_HEADSHAPE, FT_SOURCEPLOT, FT_SOURCEWRITE

% Copyright (C) 2003-2007, Robert Oostenveld
% Copyright (C) 2011-2014, Jan-Mathijs Schoffelen
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
ft_preamble loadvar functional anatomical
ft_preamble provenance functional anatomical
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% this is not supported any more as of 26/10/2011
if ischar(anatomical),
  error('please use cfg.inputfile instead of specifying the input variable as a sting');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'unused',     {'keepinside' 'voxelcoord'});
cfg = ft_checkconfig(cfg, 'deprecated', {'sourceunits', 'mriunits'});
cfg = ft_checkconfig(cfg, 'required',   'parameter');
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter', 'avg.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter', 'avg.coh', 'coh'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter', 'avg.mom', 'mom'});

% set the defaults
cfg.downsample   = ft_getopt(cfg, 'downsample', 1);
cfg.feedback     = ft_getopt(cfg, 'feedback',   'text');
cfg.interpmethod = ft_getopt(cfg, 'interpmethod', []);   % cfg.interpmethod depends on how the interpolation should be done and actual defaults will be specified below

% replace pnt by pos
anatomical = fixpos(anatomical);
functional = fixpos(functional);

% ensure the functional data to be in double precision
functional = ft_struct2double(functional);

if strcmp(cfg.interpmethod, 'nearest') && (ft_datatype(functional, 'volume+label') || ft_datatype(functional, 'source+label'))
  % the first input argument describes a parcellation or segmentation with tissue labels
  isAtlasFun = true;
else
  isAtlasFun = false;
end

if isfield(anatomical, 'transform') && isfield(anatomical, 'dim')
  % anatomical volume
  isUnstructuredAna  = false;
elseif isfield(anatomical, 'pos') && isfield(anatomical, 'dim')
  % positions that can be mapped onto a 3D regular grid
  isUnstructuredAna  = false;
elseif isfield(anatomical, 'pos')
  % anatomical data that consists of a mesh, but no smudging possible
  isUnstructuredAna  = true;
end

if isfield(functional, 'transform') && isfield(functional, 'dim')
  % functional volume
  isUnstructuredFun  = false;
elseif isfield(functional, 'pos') && isfield(functional, 'dim')
  % positions that can be mapped onto a 3D regular grid
  isUnstructuredFun  = false;
else
  isUnstructuredFun  = true;
end

if isUnstructuredAna
  anatomical = ft_checkdata(anatomical, 'datatype', {'source', 'source+label', 'mesh'}, 'inside', 'logical', 'feedback', 'yes', 'hasunit', 'yes');
else
  anatomical = ft_checkdata(anatomical, 'datatype', {'volume', 'volume+label'}, 'inside', 'logical', 'feedback', 'yes', 'hasunit', 'yes');
end

if isUnstructuredFun
  functional = ft_checkdata(functional, 'datatype', 'source', 'inside', 'logical', 'feedback', 'yes', 'hasunit', 'yes');
else
  functional = ft_checkdata(functional, 'datatype', 'volume', 'inside', 'logical', 'feedback', 'yes', 'hasunit', 'yes');
end

if ~isa(cfg.parameter, 'cell')
  cfg.parameter = {cfg.parameter};
end

% try to select all relevant parameters present in the data
if any(strcmp(cfg.parameter, 'all'))
  cfg.parameter = parameterselection('all', functional);
  for k = numel(cfg.parameter):-1:1
    % check whether the field is numeric
    tmp = getsubfield(functional, cfg.parameter{k});
    if iscell(tmp)
      cfg.parameter(k) = [];
    elseif strcmp(cfg.parameter{k}, 'pos')
      cfg.parameter(k) = [];
    end
  end
end

% ensure that the functional data has the same unit as the anatomical data
functional = ft_convert_units(functional, anatomical.unit);

if isfield(functional, 'coordsys') && isfield(anatomical, 'coordsys') && ~isequal(functional.coordsys, anatomical.coordsys)
  % FIXME is this different when smudged or not?
  % warning('the coordinate systems are not aligned');
  % error('the coordinate systems are not aligned');
end

if ~isUnstructuredAna && cfg.downsample~=1
  % downsample the anatomical volume
  tmpcfg = keepfields(cfg, {'downsample', 'showcallinfo'});
  orgcfg.parameter = cfg.parameter;
  tmpcfg.parameter = 'anatomy';
  anatomical = ft_volumedownsample(tmpcfg, anatomical);
  % restore the provenance information
  [cfg, anatomical] = rollback_provenance(cfg, anatomical);
  % restore the original parameter, it should not be 'anatomy'
  cfg.parameter = orgcfg.parameter;
end

% collect the functional volumes that should be converted
dat_name = {};
dat_array = {};
for i=1:length(cfg.parameter)
  if ~iscell(getsubfield(functional, cfg.parameter{i}))
    dat_name{end+1} = cfg.parameter{i};
    dat_array{end+1} = getsubfield(functional, cfg.parameter{i});
  else
    fprintf('not interpolating %s, since it is represented in a cell-array\n', cfg.parameter{i});
  end
end

% hmmmm, if the input data contains a time and/or freq dimension, then the output
% may be terribly blown up; most convenient would be to output only the
% smudging matrix, and project the data when plotting

if isUnstructuredFun && isUnstructuredAna && isfield(anatomical, 'orig') && isfield(anatomical.orig, 'pos') && isfield(anatomical.orig, 'tri')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % functional data defined on subset of vertices in an anatomical mesh
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % FIXME this should not be decided on the basis of the data structures but on the basis of the cfg.interpmethod option
  % FIXME the distribution of 3 geometries over the 2 structures is weird
  % FIXME a (perhaps extreme) application of this would be to interpolate data from parcels on the sheet, i.e. an inverse parcellation

  % anatomical data consists of a decimated triangulated mesh, containing
  % the original description, allowing for smudging.

  % smudge the low resolution functional data according to the strategy in
  % MNE-suite (chapter 8.3 of the manual)

  interpmat = interp_ungridded(anatomical.pos, anatomical.orig.pos, 'projmethod', 'smudge', 'triout', anatomical.orig.tri);
  interpmat(~anatomical.inside(:), :) = 0;

  % start with an empty structure, keep only some fields
  interp = keepfields(functional, {'time', 'freq'});
  interp = copyfields(anatomical, interp, {'coordsys', 'unit'});
  interp = copyfields(anatomical.orig, interp, {'pos', 'tri'});

  % identify the inside voxels after interpolation
  nzeros     = sum(interpmat~=0,2);
  newinside  = (nzeros~=0);
  newoutside = (nzeros==0);

  interp.inside = false(size(anatomical.pos,1),1);
  interp.inside(newinside) = true;

  % interpolate all functional data
  for i=1:length(dat_name)
    fprintf('interpolating %s\n', dat_name{i});

    dimord = getdimord(functional, dat_name{i});
    dimtok = tokenize(dimord, '_');
    dimf   = getdimsiz(functional, dat_name{i});
    dimf(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions

    % should be 3-D array, can have trailing singleton dimensions
    if numel(dimf)<2
      dimf(2) = 1;
    end
    if numel(dimf)<3
      dimf(3) = 1;
    end

    allav = zeros([size(anatomical.orig.pos,1), dimf(2:end)]);
    for k=1:dimf(2)
      for m=1:dimf(3)
        fv     = dat_array{i}(:,k,m);
        av     = interpmat*fv;
        av(newoutside) = nan;
        allav(:,k,m)   = av;
      end
    end
    interp = setsubfield(interp, dat_name{i}, allav);
  end


elseif isUnstructuredFun && isUnstructuredAna
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % functional data defined on a point cloud/mesh, anatomy on a volume
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set default interpmethod for this situation
  cfg.interpmethod = ft_getopt(cfg, 'interpmethod', 'nearest');
  cfg.sphereradius = ft_getopt(cfg, 'sphereradius', 0.5);
  cfg.power        = ft_getopt(cfg, 'power',        1);

  interpmat = interp_ungridded(functional.pos, anatomical.pos, 'projmethod', cfg.interpmethod, 'sphereradius', cfg.sphereradius, 'power', cfg.power); % FIXME include other key-value pairs as well
  interpmat(~anatomical.inside(:), :) = 0;

  % start with an empty structure, keep only some fields
  interp = keepfields(functional, {'time', 'freq'});
  interp = copyfields(anatomical, interp, {'pos', 'tri', 'dim', 'transform', 'coordsys', 'unit'});

  % identify the inside voxels after interpolation
  nzeros     = sum(interpmat~=0,2);
  newinside  = (nzeros~=0);
  newoutside = (nzeros==0);

  interp.inside = false(size(anatomical.pos,1),1);
  interp.inside(newinside) = true;

  % interpolate all functional data
  for i=1:length(dat_name)
    fprintf('interpolating %s\n', dat_name{i});

    dimord = getdimord(functional, dat_name{i});
    dimtok = tokenize(dimord, '_');
    dimf   = getdimsiz(functional, dat_name{i});
    dimf(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions

    % should be 3-D array, can have trailing singleton dimensions
    if numel(dimf)<2
      dimf(2) = 1;
    end
    if numel(dimf)<3
      dimf(3) = 1;
    end

    allav = zeros([size(anatomical.pos,1), dimf(2:end)]);
    for k=1:dimf(2)
      for m=1:dimf(3)
        fv     = dat_array{i}(:,k,m);
        av     = interpmat*fv;
        av(newoutside) = nan;
        allav(:,k,m)   = av;
      end
    end
    interp = setsubfield(interp, dat_name{i}, allav);
  end


elseif isUnstructuredFun && ~isUnstructuredAna
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % functional data defined on a point cloud/mesh, anatomy on a point cloud/mesh
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set default interpmethod for this situation
  cfg.interpmethod = ft_getopt(cfg, 'interpmethod', 'nearest');
  cfg.sphereradius = ft_getopt(cfg, 'sphereradius', 0.5);
  cfg.power        = ft_getopt(cfg, 'power',        1);

  [ax, ay, az] = voxelcoords(anatomical.dim, anatomical.transform);
  anatomical.pos = [ax(:) ay(:) az(:)];
  clear ax ay az

  interpmat = interp_ungridded(functional.pos, anatomical.pos, 'projmethod', cfg.interpmethod, 'sphereradius', cfg.sphereradius, 'power', cfg.power); % FIXME include other key-value pairs as well
  interpmat(~anatomical.inside(:), :) = 0;

  % start with an empty structure, keep only some fields
  interp = keepfields(functional, {'time', 'freq'});
  interp = copyfields(anatomical, interp, {'pos', 'tri', 'dim', 'transform', 'coordsys', 'unit', 'anatomy'});

  % identify the inside voxels after interpolation
  nzeros     = sum(interpmat~=0,2);
  newinside  = (nzeros~=0);
  newoutside = (nzeros==0);

  interp.inside = false(anatomical.dim);
  interp.inside(newinside) = true;

  % interpolate all functional data
  for i=1:length(dat_name)
    fprintf('interpolating %s\n', dat_name{i});

    dimord = getdimord(functional, dat_name{i});
    dimtok = tokenize(dimord, '_');
    dimf   = getdimsiz(functional, dat_name{i});
    dimf(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions

    % should be 3-D array, can have trailing singleton dimensions
    if numel(dimf)<2
      dimf(2) = 1;
    end
    if numel(dimf)<3
      dimf(3) = 1;
    end

    av    = zeros([anatomical.dim            ]);
    allav = zeros([anatomical.dim dimf(2:end)]);

    for k=1:dimf(2)
      for m=1:dimf(3)
        fv     = dat_array{i}(:,k,m);
        av(:)  = interpmat*fv;
        av(newoutside)   = nan;
        allav(:,:,:,k,m) = av;
      end
    end
    if isfield(interp, 'freq') || isfield(interp, 'time')
      % the output should be a source representation, not a volume
      allav = reshape(allav, prod(anatomical.dim), dimf(2), dimf(3));
    end
    interp = setsubfield(interp, dat_name{i}, allav);
  end

  if ~isfield(interp, 'freq') && ~isfield(interp, 'time')
      % the output should be a volume representation, not a source
    interp = rmfield(interp, 'pos');
  end

elseif ~isUnstructuredFun && isUnstructuredAna
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % functional data defined on a volume, anatomy on a point cloud/mesh
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set default interpmethod for this situation
  cfg.interpmethod = ft_getopt(cfg, 'interpmethod', 'nearest');
  cfg.sphereradius = ft_getopt(cfg, 'sphereradius', []);
  cfg.power        = ft_getopt(cfg, 'power',        1);

  % interpolate the 3D volume onto the anatomy
  if ~strcmp(cfg.interpmethod, 'project')
    % use interp_gridded
    [interpmat, dummy] = interp_gridded(functional.transform, zeros(functional.dim), anatomical.pos, 'projmethod', cfg.interpmethod, 'sphereradius', cfg.sphereradius, 'inside', functional.inside, 'power', cfg.power);

    % use interp_ungridded
    % interpmat = interp_ungridded(functional.pos, anatomical.pos, 'projmethod', cfg.interpmethod, 'sphereradius', cfg.sphereradius, 'inside', functional.inside, 'power', cfg.power);
  else
    % do the interpolation below, the current implementation of the
    % 'project' method does not output an interpmat (and is therefore quite
    % inefficient

    % set the defaults
    cfg.projvec        = ft_getopt(cfg, 'projvec',       1);
    cfg.projweight     = ft_getopt(cfg, 'projweight',    ones(size(cfg.projvec)));
    cfg.projcomb       = ft_getopt(cfg, 'projcomb',      'mean'); % or max
    cfg.projthresh     = ft_getopt(cfg, 'projthresh',    []);
  end

  % start with an empty structure, keep some fields
  interp = keepfields(functional, {'time', 'freq'});
  interp = copyfields(anatomical, interp, {'pos', 'tri', 'dim', 'transform', 'coordsys', 'unit'});

  % identify the inside voxels after interpolation
  interp.inside    = true(size(anatomical.pos,1),1);

  % interpolate all functional data
  for i=1:length(dat_name)
    fprintf('interpolating %s\n', dat_name{i});

    dimord = getdimord(functional, dat_name{i});
    dimtok = tokenize(dimord, '_');
    dimf   = getdimsiz(functional, dat_name{i});
    dimf(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions

    if prod(functional.dim)==dimf(1)
      % convert into 3-D, 4-D or 5-D array
      dimf        = [functional.dim dimf(2:end)];
      dat_array{i} = reshape(dat_array{i}, dimf);
    end

    % should be 5-D array, can have trailing singleton dimensions
    if numel(dimf)<4
      dimf(4) = 1;
    end
    if numel(dimf)<5
      dimf(5) = 1;
    end

    allav = zeros([size(anatomical.pos,1), dimf(4:end)]);
    if ~strcmp(cfg.interpmethod, 'project')
      for k=1:dimf(4)
        for m=1:dimf(5)
          fv    = dat_array{i}(:,:,:,k,m);
          fv    = fv(functional.inside(:));
          av    = interpmat*fv;
          allav(:,k,m) = av;
        end
      end
    else
      for k=1:dimf(4)
        for m=1:dimf(5)
          fv   = dat_array{i}(:,:,:,k,m);
          av   = interp_gridded(functional.transform, fv, anatomical.pos, 'dim', functional.dim, 'projmethod', 'project', 'projvec', cfg.projvec, 'projweight', cfg.projweight, 'projcomb', cfg.projcomb, 'projthresh', cfg.projthresh);
          allav(:,k,m) = av;
        end
      end
    end
    interp = setsubfield(interp, dat_name{i}, allav);
  end

elseif ~isUnstructuredFun && ~isUnstructuredAna
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % functional data defined on a volume, anatomy on a differently sampled volume
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set default interpmethod for this situation
  cfg.interpmethod = ft_getopt(cfg, 'interpmethod', 'linear');

  % start with an empty structure, keep some fields
  interp = keepfields(functional, {'time', 'freq'});
  interp = copyfields(anatomical, interp, {'pos', 'tri', 'dim', 'transform', 'coordsys', 'unit', 'anatomy'});

  % convert the anatomical voxel positions into voxel indices into the functional volume
  anatomical.transform = functional.transform \ anatomical.transform;
  functional.transform = eye(4);

  [fx, fy, fz] = voxelcoords(functional.dim, functional.transform);
  [ax, ay, az] = voxelcoords(anatomical.dim, anatomical.transform);

  % estimate the subvolume of the anatomy that is spanned by the functional volume
  minfx = 1;
  minfy = 1;
  minfz = 1;
  maxfx = functional.dim(1);
  maxfy = functional.dim(2);
  maxfz = functional.dim(3);
  sel = ax(:)>=minfx & ...
    ax(:)<=maxfx & ...
    ay(:)>=minfy & ...
    ay(:)<=maxfy & ...
    az(:)>=minfz & ...
    az(:)<=maxfz;
  fprintf('selecting subvolume of %.1f%%\n', 100*sum(sel)./prod(anatomical.dim));

  if all(functional.inside(:))
    % keep all voxels marked as inside
    interp.inside = true(anatomical.dim);
  else
    % reslice and interpolate inside
    interp.inside = zeros(anatomical.dim);
    % interpolate with method nearest
    interp.inside( sel) = my_interpn(double(functional.inside), ax(sel), ay(sel), az(sel), 'nearest', cfg.feedback);
    interp.inside(~sel) = 0;
    interp.inside = logical(interp.inside);
  end

  % prepare the grid that is used in the interpolation
  fg = [fx(:) fy(:) fz(:)];
  clear fx fy fz

  % reslice and interpolate all functional volumes
  for i=1:length(dat_name)
    fprintf('reslicing and interpolating %s\n', dat_name{i});

    dimord = getdimord(functional, dat_name{i});
    dimtok = tokenize(dimord, '_');
    dimf   = getdimsiz(functional, dat_name{i});
    dimf(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions

    if prod(functional.dim)==dimf(1)
      % convert into 3-D, 4-D or 5-D array
      dimf = [functional.dim dimf(2:end)];
      dat_array{i} = reshape(dat_array{i}, dimf);
    end

    % should be 5-D array, can have trailing singleton dimensions
    if numel(dimf)<4
      dimf(4) = 1;
    end
    if numel(dimf)<5
      dimf(5) = 1;
    end

    av    = zeros([anatomical.dim            ]);
    allav = zeros([anatomical.dim dimf(4:end)]);
    functional.inside = functional.inside(:,:,:,1,1);

    if any(dimf(4:end)>1) && ~strcmp(cfg.feedback, 'none')
      % this is needed to prevent feedback to be displayed for every time-frequency point
      warning('disabling feedback');
      cfg.feedback = 'none';
    end

    for k=1:dimf(4)
      for m=1:dimf(5)
        fv = dat_array{i}(:,:,:,k,m);
        if ~isa(fv, 'double')
          % only convert if needed, this saves memory
          fv = double(fv);
        end
        % av( sel) = my_interpn(fx, fy, fz, fv, ax(sel), ay(sel), az(sel), cfg.interpmethod, cfg.feedback);
        if islogical(dat_array{i})
          % interpolate always with method nearest
          av( sel) = my_interpn(fv, ax(sel), ay(sel), az(sel), 'nearest', cfg.feedback);
          av = logical(av);
        else
          if ~all(functional.inside(:))
            % extrapolate the outside of the functional volumes for better interpolation at the edges
            fv(~functional.inside) = griddatan(fg(functional.inside(:), :), fv(functional.inside(:)), fg(~functional.inside(:), :), 'nearest');
          end
          % interpolate functional onto anatomical grid
          av( sel) = my_interpn(fv, ax(sel), ay(sel), az(sel), cfg.interpmethod, cfg.feedback);
          av(~sel) = nan;
          av(~interp.inside) = nan;
        end
        allav(:,:,:,k,m) = av;
      end
    end
    if isfield(interp, 'freq') || isfield(interp, 'time')
      % the output should be a source representation, not a volume
      allav = reshape(allav, prod(anatomical.dim), dimf(4), dimf(5));
    end
    interp = setsubfield(interp, dat_name{i}, allav);
  end

end % computing the interpolation according to the input data

if isfield(interp, 'freq') || isfield(interp, 'time')
  % the output should be a source representation, not a volumetric representation
  if ~isfield(interp, 'pos')
    [x, y, z] = voxelcoords(interp.dim, interp.transform);
    interp.pos = [x(:) y(:) z(:)];
  end
end

if isAtlasFun
  for i=1:numel(dat_name)
    % keep the labels that describe the different tissue types
    interp = copyfields(functional, interp, [dat_name{i} 'label']);
    % replace NaNs that fall outside the labeled area with zero
    tmp = interp.(dat_name{i});
    tmp(isnan(tmp)) = 0;
    interp.(dat_name{i}) = tmp;
  end
  % remove the inside field if present
  interp = removefields(interp, 'inside');
end

if exist('interpmat', 'var')
  cfg.interpmat = interpmat;
  cfg.interpmat; % access it once to fool the cfg-tracking
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   functional anatomical
ft_postamble provenance interp
ft_postamble history    interp
ft_postamble savevar    interp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this function computes the location of all voxels in head
% coordinates in a memory efficient manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, z] = voxelcoords(dim, transform)

xgrid     = 1:dim(1);
ygrid     = 1:dim(2);
zgrid     = 1:dim(3);
npix      = prod(dim(1:2));  % number of voxels in a single slice

x = zeros(dim);
y = zeros(dim);
z = zeros(dim);
X = zeros(1,npix);
Y = zeros(1,npix);
Z = zeros(1,npix);
E = ones(1,npix);
% determine the voxel locations per slice
for i=1:dim(3)
  [X(:), Y(:), Z(:)] = ndgrid(xgrid, ygrid, zgrid(i));
  tmp = transform*[X; Y; Z; E];
  x((1:npix)+(i-1)*npix) = tmp(1,:);
  y((1:npix)+(i-1)*npix) = tmp(2,:);
  z((1:npix)+(i-1)*npix) = tmp(3,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for memory efficient interpolation
% the only reason for this function is that it does the interpolation in smaller chuncks
% this prevents memory problems that I often encountered here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [av] = my_interpn(fx, fy, fz, fv, ax, ay, az, interpmethod, feedback);
function [av] = my_interpn(fv, ax, ay, az, interpmethod, feedback)

num = numel(ax);            % total number of voxels
blocksize = floor(num/20);  % number of voxels to interpolate at once, split it into 20 chuncks
lastblock = 0;              % boolean flag for while loop
sel = 1:blocksize;          % selection of voxels that are interpolated, this is the first chunck
av  = zeros(size(ax));
ft_progress('init', feedback, 'interpolating');
while (1)
  ft_progress(sel(1)/num, 'interpolating %.1f%%\n', 100*sel(1)/num);
  if sel(end)>=num
    sel = sel(1):num;
    lastblock = 1;
  end
  av(sel) = interpn(fv, ax(sel), ay(sel), az(sel), interpmethod);
  if lastblock
    break
  end
  sel = sel + blocksize;
end
ft_progress('close');
