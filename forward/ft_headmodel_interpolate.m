function headmodel = ft_headmodel_interpolate(filename, sens, sourcemodel, varargin)

% FT_HEADMODEL_INTERPOLATE describes a volume conduction model of the head in which
% subsequent leadfield computations can be performed using a simple interpolation
% scheme.
%
% Use as
%   headmodel = ft_headmodel_interpolate(filename, sens, leadfield)
% or
%   headmodel = ft_headmodel_interpolate(filename, sens, leadfield)
%
% The input parameters are the filename to which the model will be written,
% the electrode definition (see ft_DATATYPE_SENS). The third input argument
% is either a pre-computed leadfield structure from FT_PREPARE_LEADFIELD
% or a the output of a previous call to FT_HEADMODEL_INTERPOLATE.
%
% The output volume conduction model is stored on disk in a MATLAB file together with a
% number of NIFTi files. The mat file contains a structure with the following fields
%   headmodel.sens        = structure, electrode sensor description, see FT_DATATYE_SENS
%   headmodel.filename    = cell-array with NIFTI filenames, one file per channel
% and contains
%   headmodel.dim         = [Nx Ny Nz] vector with the number of grid points along each dimension
%   headmodel.transform   = 4x4 homogenous transformation matrix
%   headmodel.unit        = string with the geometrical units of the positions, e.g. 'cm' or 'mm'
% to describe the source positions.
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if ischar(sourcemodel)
  % the input is a filename that points to a BESA precomputed leadfield
  filename = sourcemodel;
  clear sourcemodel
  
  % this requires the BESA functions
  ft_hastoolbox('besa', 1);
  
  % get the filename with the path but without the extension
  [p, f, x] = fileparts(filename);
  if isempty(p)
    p = pwd;
  end
  lftfile = fullfile(p, [f, '.lft']);
  locfile = fullfile(p, [f, '.loc']);
  
  % Read source positions
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
  
  insideindx = sub2ind(dim, ind(1,:), ind(2,:), ind(3,:)); % note that ind is transposed
  
  if false
    % this shows how the positions are reconstructed from dim+transform+insideindx
    [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
    vox = [X(:) Y(:) Z(:)];
    head = ft_warp_apply(transform, vox);
    assert(norm(head(insideindx,:)-ssg)/norm(ssg)<1e-9); % there is a little bit rounding off error
  end
  
  sourcemodel           = [];
  sourcemodel.dim       = dim;
  sourcemodel.transform = transform;
  sourcemodel.inside    = false(prod(dim),1);
  sourcemodel.inside(insideindx) = true;
  sourcemodel.leadfield = cell(dim);
  
  % ensure that it has geometrical units (probably mm)
  sourcemodel = ft_determine_units(sourcemodel);
  
  % Read leadfield, all channels, all locations, 3 orientations
  [lftdim, lft] = readBESAlft(lftfile);
  
  assert(lftdim(1)==length(sens.label), 'inconsistent number of electrodes');
  assert(lftdim(2)==length(insideindx), 'inconsistent number of source positions');
  assert(lftdim(3)==3, 'unexpected number of leadfield columns');
  assert(isequal(sourcemodel.unit, sens.unit), 'inconsistent geometrical units');
  
  
  for i=1:length(insideindx)
    sel = 3*(i-1)+(1:3);
    sourcemodel.leadfield{insideindx(i)} = lft(:,sel);
  end
  
  fprintf('finished import of BESA leadfield file\n');
end % process the BESA file, sourcemodel is now compatible with FT_PREPARE_LEADFIELD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART TWO: write the leadfield to a set of nifti files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(sourcemodel, 'leadfield')
  % the input pre-computed leadfields reflect the output of FT_PREPARE_LEADFIELD
  % which should be reorganized into channel-specific volumes and stored to disk as nifti files
  
  % ensure that it is represented as 3-D volume
  sourcemodel = ft_checkdata(sourcemodel, 'datatype', 'volume');
  
  nchan = length(sens.label);
  if size(sourcemodel.leadfield{insideindx(1)},1)~=nchan
    ft_error('the number of channels does not match');
  end
  
  headmodel = [];
  headmodel.type      = 'interpolate';
  headmodel.dim       = sourcemodel.dim;
  headmodel.transform = sourcemodel.transform;
  headmodel.inside    = false(sourcemodel.dim);
  headmodel.inside(insideindx) = true;
  headmodel.sens      = sens;
  headmodel.filename  = cell(size(sens.label));
  
  if isfield(sourcemodel, 'unit')
    % get the units from the sourcemodel
    headmodel.unit = sourcemodel.unit;
  else
    % estimate the units
    headmodel = ft_determine_units(headmodel);
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

  % these indices only have to be determined once and speed-up the reassignment
  [ind1, ind2, ind3] = ndgrid(1:headmodel.dim(1), 1:headmodel.dim(2), 1:headmodel.dim(3));
  
  for i=1:nchan
    dat = zeros([headmodel.dim 3]);
    for j=insideindx(:)'
      % [i1, i2, i3] = ind2sub(headmodel.dim, j);
      % ind2sub is slow, simply look them up instead
      i1 = ind1(j);
      i2 = ind2(j);
      i3 = ind3(j);
      dat(i1, i2, i3, :) = sourcemodel.leadfield{j}(i,:);
    end
    
    if istrue(smooth)
      if i == 1
        ft_write_mri(masklf,~~dat(:, :, :, 1), 'transform', sourcemodel.transform, 'spmversion', 'SPM12', 'dataformat', 'nifti_spm');
        spm_smooth(masklf, smasklf, sourcemodel.transform(1,1)*[1 1 1]);
        mask = spm_read_vols(spm_vol(smasklf));
        
        spm_unlink(masklf);
        spm_unlink(smasklf);
      end
      
      ft_write_mri(rawlf, dat, 'transform', sourcemodel.transform, 'spmversion', 'SPM12', 'dataformat', 'nifti_spm');
      spm_smooth(rawlf, srawlf, sourcemodel.transform(1,1)*[1 1 1]);
      dat = spm_read_vols(spm_vol(srawlf));
      dat = dat./repmat(mask, [1 1 1, size(dat, 4)]);
      dat(~isfinite(dat)) = 0;
      
      spm_unlink(rawlf);
      spm_unlink(srawlf);
    end
    
    
    headmodel.filename{i} = sprintf('%s_%s.nii', filename, sens.label{i});
    fprintf('writing single channel leadfield to %s\n', headmodel.filename{i})
    
    if exist('spm_bsplinc', 'file')
      dat = cat(4, dat, 0*dat);
      for k = 1:3
        dat(:, :, :, k+3) = spm_bsplinc(dat(:, :, :, k), [4 4 4 0 0 0]);
      end
    end
    
    ft_write_mri(headmodel.filename{i}, dat , 'transform', sourcemodel.transform, 'spmversion', 'SPM12', 'dataformat', 'nifti_spm');
    
  end
  
  filename = sprintf('%s.mat', filename);
  fprintf('writing volume conduction model metadata to %s\n', filename)
  save(filename, 'headmodel');
  
elseif isfield(sourcemodel, 'filename')
  % the input pre-computed leadfields reflect the output of FT_HEADMODEL_INTERPOLATE,
  % which should be re-interpolated on the channel level and then stored to disk as nifti files
  ft_hastoolbox('spm8up', 1);
  
  inputvol = sourcemodel;
  
  if ~isfield(sens, 'tra')
    sens.tra = eye(length(sens.label));
  end
  
  if ~isfield(inputvol.sens, 'tra')
    inputvol.sens.tra = eye(length(inputvol.sens.label));
  end
  
  % create a 2D projection and triangulation
  pos = inputvol.sens.elecpos;
  prj = elproj(pos);
  tri = delaunay(prj(:,1), prj(:,2));
  
  % project the electrodes on the triangulation and compute the
  % bilinear interpolation from the original to the new electrodes
  [el, prj] = project_elec(sens.elecpos, pos, tri);
  tra = transfer_elec(pos, tri, el);
  
  % define the spaces and the number of elements that they comprise
  n1 = length(inputvol.sens.label);    % computed channels
  n2 = size(inputvol.sens.elecpos,1);  % computed electrode positions
  n3 = size(sens.elecpos,1);      % desired electrode positions
  n4 = length(sens.label);        % desired channels
  
  % this is the montage for getting the the desired channels from the desired electrode positions
  make4from3.labelold = cell(n3,1);
  make4from3.labelnew = sens.label;
  make4from3.tra      = sens.tra;
  for i=1:n3
    make4from3.labelold{i} = sprintf('3to4_%d', i);
  end
  
  % this is the montage for getting the computed channels from the computed electrode positions
  make1from2.labelold = cell(n2,1);
  make1from2.labelnew = inputvol.sens.label;
  make1from2.tra      = inputvol.sens.tra;
  for i=1:n2
    make1from2.labelold{i} = sprintf('2to1_%d', i);
  end
  
  % this is the montage that maps the computed electrode positions to the desired positions
  make3from2.labelold = make1from2.labelold; % the computed electrodes
  make3from2.labelnew = make4from3.labelold; % the desired electrodes
  make3from2.tra      = tra;
  
  % the following should be read as a sequence of left-hand multiplications
  % we need           make4from1
  % we can make it as make4from3 * make3from2 * make2from1
  % or as             make4from3 * make3from2 * inv(make1from2)
  
  make4from1.tra      = make4from3.tra * make3from2.tra / make1from2.tra;
  make4from1.labelold = make1from2.labelnew;
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
        dat(:, :, :, k+3) = spm_bsplinc(dat(:, :, :, k), [4 4 4 0 0 0]);
      end
    end
    outputvol.filename{i} = sprintf('%s_%s.nii', filename, sens.label{i});
    fprintf('writing single channel leadfield to %s\n', outputvol.filename{i})
    ft_write_mri(outputvol.filename{i}, dat, 'transform', outputvol.transform, 'spmversion', 'SPM12', 'dataformat', 'nifti_spm');
  end
  
  % update the volume conductor
  outputvol.sens = sens;
  
  % rename and save to disk
  headmodel = outputvol;
  clear inputvol outputvol
  
  filename = sprintf('%s.mat', filename);
  fprintf('writing volume conductor structure to %s\n', filename)
  save(filename, 'headmodel');
  
end


