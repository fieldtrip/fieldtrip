function vol = ft_headmodel_interpolate(filename, sens, grid, varargin)

% FT_HEADMODEL_INTERPOLATE describes a volume conduction model of the head in which
% subsequent leadfield computations can be performed using a simple interpolation
% scheme.
%
% Use as
%   vol = ft_headmodel_interpolate(filename, sens, leadfield)
% where the input parameters are the filename to which the model has to be written, sens
% contains the electrode definition (see ft_DATATYPE_SENS) and leadfield contains the
% pre-computed leadfielda (see FT_PRERARE_LEADFIELD)
%
% The output volume conduction model is stored on disk in a MATLAB file together with a
% number of NIFTi files. The mat file contains a structure with the following fields
%   vol.sens        = structure, electrode sensor description, see FT_DATATYE_SENS
%   vol.filename    = cell-array with NIFTi filenames, one per channel
% and contains either
%   vol.dim         = [Nx Ny Nz] vector with the number of grid points along each dimension
%   vol.transform   = 4x4 homogenous transformation matrix
% or
%   vol.pos         = Nx3 matrix with source positions
% to describe the source positions.
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% check the validity of the input arguments
assert(ft_datatype(sens, 'sens'), 'the second input argument should be a sensor definition');
grid = ft_checkdata(grid, 'datatype', 'volume', 'inside', 'index');

nchan = length(sens.label);

[p, f, x] = fileparts(filename);
filename = fullfile(p, f); % without the extension

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
  lf = cat(4, lfx, lfy, lfz);
  vol.filename{i} = sprintf('%s_%s.nii', filename, sens.label{i});
  fprintf('writing single channel leadfield to %s\n', vol.filename{i})
  ft_write_mri(vol.filename{i}, lf);
end

filename = sprintf('%s.mat', filename);
fprintf('writing volume conductor structure to %s\n', filename)
save(filename, 'vol');

