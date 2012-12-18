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

% get the optional input arguments
conductivity = ft_getopt(varargin, 'conductivity'); % default is determined below
fitind       = ft_getopt(varargin, 'fitind', 'all');
unit         = ft_getopt(varargin, 'unit');

% start with an empty volume conductor
vol = [];

if ~isempty(unit)
  % use the user-specified units for the output
  vol.unit = unit;
elseif isfield(geometry, 'unit')
  % copy the geometrical units into he volume conductor
  % assume that in case of multiple meshes that they have the same units
  vol.unit = geometry(1).unit;
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


vol.type = 'concentricspheres';
vol.o    = single_o;              % specify the center of the spheres
vol.c    = conductivity;          % specify the conductivity of the spheres
vol      = ft_convert_units(vol); % ensure the object to have a unit

% sort the spheres from the smallest to the largest ('insidefirst' order)
[vol.r, indx] = sort(vol.r);

if isempty(vol.c)
  % it being empty indicates that the user did not specify a conductivity, use a default instead
  if length(vol.r)==1
    vol.c = 1; % brain
  elseif length(vol.r)==3
    vol.c = [0.3300   0.0042 0.3300]; % brain,      skull, skin
  elseif length(vol.r)==4
    vol.c = [0.3300 1 0.0042 0.3300]; % brain, csf, skull, skin
  else
    error('conductivity values should be specified for each tissue type');
  end
else
  % the conductivity as specified by the user should be in the same order as the geometries
  % sort the spheres from the smallest to the largest ('insidefirst' order)
  vol.c = vol.c(indx);
end

for i=1:numel(geometry)
  fprintf('concentric sphere %d: radius = %.1f, conductivity = %f\n', i, vol.r(i), vol.c(i));
end

