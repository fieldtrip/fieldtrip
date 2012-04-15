function vol = ft_headmodel_concentricspheres(geometry, varargin)

% FT_HEADMODEL_CONCENTRICSPHERES creates a volume conduction model
% of the head based on three or four concentric spheres. For a 3-sphere
% model the spheres represent the skin surface, the outside of the
% skull and the inside of the skull For a 4-sphere model, the surfaces
% describe the skin, the outside-skull, the inside-skull and the inside of the
% cerebro-spinal fluid (CSF) boundaries.
%
% The innermost surface is sometimes also referred to as the brain
% surface, i.e. as the outside of the brain volume.
%
% This function takes as input a single headshape described with
% points and fits the spheres to this surface. If you have a set of
% points describing each surface, then this function fits the spheres
% to all individual surfaces.
%
% Use as
%   vol = ft_headmodel_concentricspheres(geometry, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   conductivity     = vector with the conductivity of each compartment
%   fitind           = vector with indices of the surfaces to use in fitting the center of the spheres
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% get the optional input arguments
conductivity = ft_getopt(varargin, 'conductivity');
fitind       = ft_getopt(varargin, 'fitind', 'all');
unit         = ft_getopt(varargin,'unit');

% start with an empty volume conductor
vol = [];

if ~isempty(unit)
  % use the user-specified units for the output
  vol.unit = geometry.unit;
elseif isfield(geometry, 'unit')
  % copy the geometrical units into he volume conductor
  vol.unit = geometry.unit;
end

if isnumeric(geometry) && size(geometry,2)==3
  % assume that it is a Nx3 array with vertices
  geometry.pnt = geometry;
elseif isstruct(geometry) && isfield(geometry,'bnd')
  % take the triangulated surface
  geometry = geometry.bnd;
end

if isequal(fitind, 'all')
  fitind = 1:numel(geometry);
end

% determine the number of compartments
numboundaries = numel(geometry);

if isempty(conductivity)
  warning('No conductivity is declared, Assuming standard values\n')
  if numboundaries == 1
    conductivity = 1;
  elseif numboundaries == 3
    % skin/skull/brain
    conductivity = [1 1/80 1] * 0.33;
  elseif numboundaries == 4
    %FIXME: check for better default values here
    % skin / outer skull / inner skull / brain
    conductivity = [1 1/80 1 1] * 0.33;
  else
    error('Conductivity values are required!')
  end
end

if numel(conductivity)~=numboundaries
  error('a conductivity value should be specified for each compartment');
else
  % assign the conductivity of each compartment
  vol.cond = conductivity;
end

% concatenate the vertices of all surfaces
pnt = [];
for i = fitind
  pnt = [pnt ; geometry(i).pnt];
end

% remove double vertices
pnt  = unique(pnt, 'rows');
npnt = size(pnt, 1);

% fit a single sphere to all combined headshape points
[single_o, single_r] = fitsphere(pnt);
fprintf('initial sphere: number of unique surface points = %d\n', npnt);
fprintf('initial sphere: center = [%.1f %.1f %.1f]\n', single_o(1), single_o(2), single_o(3));
fprintf('initial sphere: radius = %.1f\n', single_r);

% fit the radius of each concentric sphere to the corresponding surface points
for i = 1:numel(geometry)
  npnt     = size(geometry(i).pnt,1);
  dist     = sqrt(sum(((geometry(i).pnt - repmat(single_o, npnt, 1)).^2), 2));
  vol.r(i) = mean(dist);
end

% specify the center of the spheres
vol.o    = single_o;
vol.c    = conductivity;
vol.type = 'concentric';
vol      = ft_convert_units(vol); % ensure the object to have a unit

% sort the spheres from the smallest to the largest ('insidefirst' order)
[vol.r, indx] = sort(vol.r);
vol.c = vol.c(indx);

for i=1:numel(geometry)
  fprintf('concentric sphere %d: radius = %.1f, conductivity = %f\n', i, vol.r(i), vol.c(i));
end

