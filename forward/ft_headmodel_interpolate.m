function vol = ft_headmodel_interpolate(filename, sens, grid, varargin)

% FT_HEADMODEL_INTERPOLATE describes a volume conduction model of the head in which
% subsequent leadfield computations can be performed using a simple interpolation
% scheme.
%
% Use as
%   vol = ft_headmodel_interpolate(filename, sens, leadfield)
% or
%   vol = ft_headmodel_interpolate(filename, sens, leadfield)
%
% The input parameters are the filename to which the model will be written,
% the electrode definition (see ft_DATATYPE_SENS). The third input argument
% is either a pre-computed leadfield structure from FT_PREPARE_LEADFIELD
% or a the output of a previous call to FT_HEADMODEL_INTERPOLATE.
%
% The output volume conduction model is stored on disk in a MATLAB file together with a
% number of NIFTi files. The mat file contains a structure with the following fields
%   vol.sens        = structure, electrode sensor description, see FT_DATATYE_SENS
%   vol.filename    = cell-array with NIFTI filenames, one file per channel
% and contains
%   vol.dim         = [Nx Ny Nz] vector with the number of grid points along each dimension
%   vol.transform   = 4x4 homogenous transformation matrix
%   vol.unit        = string with the geometrical units of the positions, e.g. 'cm' or 'mm'
% to describe the source positions.
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

% check the validity of the input arguments
assert(ft_datatype(sens, 'sens'), 'the second input argument should be a sensor definition');

% get the optional input arguments
smooth = ft_getopt(varargin, 'smooth', true);

% get the filename with the path but without the extension
[p, f, x] = fileparts(filename);
if isempty(p)
  p = pwd;
end
filename = fullfile(p, f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART ONE (optional), read the pre-computed besa leadfield
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(grid)
  % the input is a filename that points to a BESA precomputed leadfield
  hdmfile = grid;
  clear grid
  
  % this requires the BESA functions
  ft_hastoolbox('besa', 1);
  
  % get the filename with the path but without the extension
  [p, f, x] = fileparts(hdmfile);
  if isempty(p)
    p = pwd;
  end
  lftfile = fullfile(p, [f, '.lft']);
  locfile = fullfile(p, [f, '.loc']);
  
  % Read source space grid nodes
  [ssg, IdxNeighbour] = readBESAloc(locfile);
  fprintf('Number of nodes: %i\n', size(ssg, 1));
  fprintf('Number of neighbours/node: %i\n', size(IdxNeighbour, 1));
  
  % the locations are represented as a Nx3 list of grid points
  % convert the representation to dim/transform/inside
  
  minx = min(ssg(:,1));
  maxx = max(ssg(:,1));
  miny = min(ssg(:,2));
  maxy = max(ssg(:,2));
  minz = min(ssg(:,3));
  maxz = max(ssg(:,3));
  
  xgrid = sort(unique(ssg(:,1)));
  ygrid = sort(unique(ssg(:,2)));
  zgrid = sort(unique(ssg(:,3)));
  
  dim = [length(xgrid) length(ygrid) length(zgrid)];
  
  pos = ssg;
  ind = zeros(size(pos));
  
  for i=1:size(ssg,1)
    ind(i,1) = find(xgrid==pos(i,1));
    ind(i,2) = find(ygrid==pos(i,2));
    ind(i,3) = find(zgrid==pos(i,3));
  end
  
  ind = ind';ind(4,:) = 1;
  pos = pos';pos(4,:) = 1;
  transform = pos/ind;
  
  inside = sub2ind(dim, ind(1,:), ind(2,:), ind(3,:)); % note that ind is transposed
  
  if false
    % this shows how the positions are reconstructed from dim+transform+inside
    [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
    vox = [X(:) Y(:) Z(:)];
    head = ft_warp_apply(transform, vox);
    assert(norm(head(inside,:)-ssg)/norm(ssg)<1e-9); % there is a little bit rounding off error
  end
  
  grid           = [];
  grid.dim       = dim;
  grid.transform = transform;
  grid.inside    = inside; % all other grid points are assumed to be "outside"
  grid.leadfield = cell(dim);
  
  % ensure that it has geometrical units (probably mm)
  grid = ft_convert_units(grid);
  
  % Read leadfield, all channels, all locations, 3 orientations
  [lftdim, lft] = readBESAlft(lftfile);
  
  assert(lftdim(1)==length(sens.label), 'inconsistent number of electrodes');
  assert(lftdim(2)==length(inside), 'inconsistent number of grid positions');
  assert(lftdim(3)==3, 'unexpected number of leadfield columns');
  assert(isequal(grid.unit, sens.unit), 'inconsistent geometrical units');
  
  
  for i=1:length(grid.inside)
    sel = 3*(i-1)+(1:3);
    grid.leadfield{grid.inside(i)} = lft(:,sel);
  end
  
  fprintf('finished import of BESA leadfield file\n');
end % process the BESA file, grid is now compatible with FT_PREPARE_LEADFIELD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART TWO: write the leadfield to a set of nifti files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(grid, 'leadfield')
  % the input pre-computed leadfields reflect the output of FT_PREPARE_LEADFIELD
  % which should be reorganized into channel-specific volumes and stored to disk as nifti files
  
  % ensure that it is represented as 3-D volume
  grid = ft_checkdata(grid, 'datatype', 'volume', 'inside', 'index');
  
  nchan = length(sens.label);
  if size(grid.leadfield{grid.inside(1)},1)~=nchan
    error('the number of channels does not match');
  end
  
  vol = [];
  vol.type      = 'interpolate';
  vol.dim       = grid.dim;
  vol.transform = grid.transform;
  vol.inside    = grid.inside;
  vol.sens      = sens;
  vol.filename  = cell(size(sens.label));
  
  if isfield(grid, 'unit')
    % get the units from the dipole grid
    vol.unit = grid.unit;
  else
    % estimate the units
    vol = ft_convert_units(vol);
  end
  
  % these go in the same directory as the other nii files, they will be removed after use
  masklf  = fullfile(p, 'masklf.nii');
  smasklf = fullfile(p, 'smasklf.nii');
  rawlf   = fullfile(p, 'rawlf.nii');
  srawlf  = fullfile(p, 'srawlf.nii');
  
  % ensure that the output directory exists
  if ~exist(p, 'dir')
    mkdir(p);
  end
  
  for i=1:nchan
    dat = zeros([vol.dim 3]);
    for j=grid.inside(:)'
      [i1, i2, i3] = ind2sub(vol.dim, j);
      dat(i1, i2, i3, 1) = grid.leadfield{j}(i,1);
      dat(i1, i2, i3, 2) = grid.leadfield{j}(i,2);
      dat(i1, i2, i3, 3) = grid.leadfield{j}(i,3);
    end
    
    if istrue(smooth)
      if i == 1
        ft_write_mri(masklf,~~dat(:, :, :, 1), 'transform', grid.transform, 'spmversion', 'SPM12', 'dataformat', 'nifti_spm');
        spm_smooth(masklf, smasklf, grid.transform(1,1)*[1 1 1]);
        mask = spm_read_vols(spm_vol(smasklf));
        
        spm_unlink(masklf);
        spm_unlink(smasklf);
      end
      
      ft_write_mri(rawlf, dat, 'transform', grid.transform, 'spmversion', 'SPM12', 'dataformat', 'nifti_spm');
      spm_smooth(rawlf, srawlf, grid.transform(1,1)*[1 1 1]);
      dat = spm_read_vols(spm_vol(srawlf));
      dat = dat./repmat(mask, [1 1 1, size(dat, 4)]);
      dat(~isfinite(dat)) = 0;
      
      spm_unlink(rawlf);
      spm_unlink(srawlf);
    end
    
    
    vol.filename{i} = sprintf('%s_%s.nii', filename, sens.label{i});
    fprintf('writing single channel leadfield to %s\n', vol.filename{i})
    
    if exist('spm_bsplinc', 'file')
      dat = cat(4, dat, 0*dat);
      for k = 1:3
        dat(:, :, :, k+3) = spm_bsplinc(squeeze(dat(:, :, :, k)), [4 4 4 0 0 0]);
      end
    end
    
    ft_write_mri(vol.filename{i}, dat , 'transform', grid.transform, 'spmversion', 'SPM12', 'dataformat', 'nifti_spm');
    
  end
  
  filename = sprintf('%s.mat', filename);
  fprintf('writing volume conduction model metadata to %s\n', filename)
  save(filename, 'vol');
  
elseif isfield(grid, 'filename')
  % the input pre-computed leadfields reflect the output of FT_HEADMODEL_INTERPOLATE,
  % which should be re-interpolated on the channel level and then stored to disk as nifti files
  ft_hastoolbox('spm8up', 1);
  
  inputvol = grid;
  
  if ~isfield(sens, 'tra')
    sens.tra = eye(length(sens.label));
  end
  
  if ~isfield(inputvol.sens, 'tra')
    inputvol.sens.tra = eye(length(inputvol.sens.label));
  end
  
  % create a 2D projection and triangulation
  pnt = inputvol.sens.elecpos;
  prj = elproj(pnt);
  tri = delaunay(prj(:,1), prj(:,2));
  
  % project the electrodes on the triangulation and compute the
  % bilinear interpolation from the original to the new electrodes
  [el, prj] = project_elec(sens.elecpos, pnt, tri);
  tra = transfer_elec(pnt, tri, el);
  
  % define the spaces and the number of elements that they comprise
  n1 = length(inputvol.sens.label);    % computed channels
  n2 = size(inputvol.sens.elecpos,1);  % computed electrode positions
  n3 = size(sens.elecpos,1);      % desired electrode positions
  n4 = length(sens.label);        % desired channels
  
  % this is the montage for getting the the desired channels from the desired electrode positions
  make4from3.labelorg = cell(n3,1);
  make4from3.labelnew = sens.label;
  make4from3.tra      = sens.tra;
  for i=1:n3
    make4from3.labelorg{i} = sprintf('3to4_%d', i);
  end
  
  % this is the montage for getting the computed channels from the computed electrode positions
  make1from2.labelorg = cell(n2,1);
  make1from2.labelnew = inputvol.sens.label;
  make1from2.tra      = inputvol.sens.tra;
  for i=1:n2
    make1from2.labelorg{i} = sprintf('2to1_%d', i);
  end
  
  % this is the montage that maps the computed electrode positions to the desired positions
  make3from2.labelorg = make1from2.labelorg; % the computed electrodes
  make3from2.labelnew = make4from3.labelorg; % the desired electrodes
  make3from2.tra      = tra;
  
  % the following should be read as a sequence of left-hand multiplications
  % we need           make4from1
  % we can make it as make4from3 * make3from2 * make2from1
  % or as             make4from3 * make3from2 * inv(make1from2)
  
  make4from1.tra      = make4from3.tra * make3from2.tra / make1from2.tra;
  make4from1.labelorg = make1from2.labelnew;
  make4from1.labelnew = make4from3.labelnew;
  
  sens = ft_apply_montage(inputvol.sens, make4from1, 'keepunused', 'no');
  % make the rounding off errors equal to zero
  sens.tra(sens.tra<10*eps) = 0;
  
  % map the leadfields for the old channels into memory
  chan = cell(1,n1);
  for j=1:n1
    chan{j} = nifti(inputvol.filename{j});
  end
  
  % construct the new volume conduction model
  outputvol.type      = inputvol.type;
  outputvol.dim       = inputvol.dim;
  outputvol.transform = inputvol.transform;
  outputvol.inside    = inputvol.inside;
  outputvol.sens      = sens;
  outputvol.filename  = cell(1,n4);
  outputvol.unit      = inputvol.unit;
  
  for i=1:n4 % each of the new channels
    dat = zeros([outputvol.dim 3]);
    for j=1:n1  % each of the old channels
      weight = make4from1.tra(i,j);
      if weight
        % interpolate the leadfields from the old to the new channels
        dat = dat + weight * chan{j}.dat(:,:,:,1:3);
      end
    end
    if exist('spm_bsplinc', 'file')
      dat = cat(4, dat, 0*dat);
      for k = 1:3
        dat(:, :, :, k+3) = spm_bsplinc(squeeze(dat(:, :, :, k)), [4 4 4 0 0 0]);
      end
    end
    outputvol.filename{i} = sprintf('%s_%s.nii', filename, sens.label{i});
    fprintf('writing single channel leadfield to %s\n', outputvol.filename{i})
    ft_write_mri(outputvol.filename{i}, dat, 'transform', outputvol.transform, 'spmversion', 'SPM12', 'dataformat', 'nifti_spm');
  end
  
  % update the volume conductor
  outputvol.sens = sens;
  
  % rename and save to disk
  vol = outputvol;
  clear inputvol outputvol
  
  filename = sprintf('%s.mat', filename);
  fprintf('writing volume conductor structure to %s\n', filename)
  save(filename, 'vol');
  
end
