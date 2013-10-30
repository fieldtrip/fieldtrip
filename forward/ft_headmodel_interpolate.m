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

% the file with the path but without the extension
[p, f, x] = fileparts(filename);

if isempty(p)
    p = pwd;
end

res = mkdir(p, f);

filename = fullfile(p, f, f);

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
  
  lfx = zeros(vol.dim);
  lfy = zeros(vol.dim);
  lfz = zeros(vol.dim);
  
  for i=1:nchan
    for j=grid.inside(:)'
      lfx(j) = grid.leadfield{j}(i,1);
      lfy(j) = grid.leadfield{j}(i,2);
      lfz(j) = grid.leadfield{j}(i,3);
    end
    dat = cat(4, lfx, lfy, lfz);
    if exist('spm_bsplinc', 'file')
        dat = cat(4, dat, 0*dat);
        for k = 1:3
            dat(:, :, :, k+3) = spm_bsplinc(squeeze(dat(:, :, :, k)), [4 4 4 0 0 0]);
        end
    end
    vol.filename{i} = sprintf('%s_%s.nii', filename, sens.label{i});
    fprintf('writing single channel leadfield to %s\n', vol.filename{i})
    ft_write_mri(vol.filename{i}, dat, 'spmversion', 'SPM12', 'dataformat', 'nifti_spm');
  end
  
  filename = sprintf('%s.mat', filename);
  fprintf('writing volume conductor structure to %s\n', filename)
  save(filename, 'vol');
  
elseif isfield(grid, 'filename')
  % the input pre-computed leadfields reflect the output of FT_HEADMODEL_INTERPOLATE,
  % which should be re-interpolated on the channel level and then stored to disk as nifti files
  ft_hastoolbox('spm8up', 1);
  
  inputvol = grid;
  
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
    make4from3.labelorg{i} = sprintf('p3_%d', i);
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
