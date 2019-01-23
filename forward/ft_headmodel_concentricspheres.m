function headmodel = ft_headmodel_concentricspheres(mesh, varargin)

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
%   headmodel = ft_headmodel_concentricspheres(mesh, ...)
%
% Optional input arguments should be specified in key-value pairs and can include
%   conductivity = vector with the conductivity of each compartment
%   fitind       = vector with indices of the surfaces to use in fitting the center of the spheres
%   order        = number of iterations in series expansion (default = 60)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2012-2018, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

% get the optional input arguments
conductivity = ft_getopt(varargin, 'conductivity'); % default is determined below
fitind       = ft_getopt(varargin, 'fitind', 'all');
order        = ft_getopt(varargin, 'order', 60);

if any(strcmp(varargin(1:2:end), 'unit')) || any(strcmp(varargin(1:2:end), 'units'))
  % the geometrical units should be specified in the input mesh
  ft_error('the ''unit'' option is not supported any more');
end

if isnumeric(mesh) && size(mesh,2)==3
  % assume that it is a Nx3 array with vertices
  % convert it to a structure, this is needed to determine the units further down
  mesh = struct('pos', mesh);
elseif isstruct(mesh) && isfield(mesh, 'bnd')
  % take the triangulated surfaces from the input structure
  mesh = mesh.bnd;
end

% replace pnt with pos
mesh = fixpos(mesh);

if ~isstruct(mesh) || ~isfield(mesh, 'pos')
  ft_error('the input mesh should be a set of points or a single triangulated surface')
end

% start with an empty volume conductor
headmodel = [];

% ensure that the mesh has units, estimate them if needed
mesh = ft_determine_units(mesh);

% copy the geometrical units into the volume conductor
headmodel.unit = mesh(1).unit;

if isequal(fitind, 'all')
  fitind = 1:numel(mesh);
end

% concatenate the vertices of all surfaces
pos = {mesh.pos};
pos = cat(1, pos{:});

% remove double vertices
pos  = unique(pos, 'rows');
npos = size(pos, 1);

% fit a single sphere to all combined headshape points
[single_o, single_r] = fitsphere(pos);
fprintf('initial sphere: number of unique surface points = %d\n', npos);
fprintf('initial sphere: center = [%.1f %.1f %.1f]\n', single_o(1), single_o(2), single_o(3));
fprintf('initial sphere: radius = %.1f\n', single_r);

% fit the radius of each concentric sphere to the corresponding surface points
for i = 1:numel(mesh)
  npos     = size(mesh(i).pos,1);
  dist     = sqrt(sum(((mesh(i).pos - repmat(single_o, npos, 1)).^2), 2));
  headmodel.r(i) = mean(dist);
end

headmodel.o    = single_o;              % specify the center of the spheres
headmodel.cond = conductivity;          % specify the conductivity of the spheres, can be empty up to here
headmodel.type = 'concentricspheres';

% sort the spheres from the smallest to the largest
[headmodel.r, indx] = sort(headmodel.r);

if isempty(headmodel.cond)
  % it being empty indicates that the user did not specify a conductivity, use a default instead
  if length(headmodel.r)==1
    headmodel.cond = 1;                        % brain
  elseif length(headmodel.r)==3
    headmodel.cond = [0.3300   0.0042 0.3300]; % brain,      skull, skin
  elseif length(headmodel.r)==4
    headmodel.cond = [0.3300 1 0.0042 0.3300]; % brain, csf, skull, skin
  else
    ft_error('conductivity values should be specified for each tissue type');
  end
else
  % the conductivity as specified by the user should be in the same order as the geometries
  % sort the spheres from the smallest to the largest ('insidefirst' order)
  headmodel.cond = headmodel.cond(indx);
end

for i=1:numel(mesh)
  fprintf('concentric sphere %d: radius = %f, conductivity = %f\n', i, headmodel.r(i), headmodel.cond(i));
end

% replicate the spheres if there are fewer than four
switch numel(mesh)
  case 1
    headmodel.r    = [headmodel.r(1)    headmodel.r(1)    headmodel.r(1)    headmodel.r(1)];
    headmodel.cond = [headmodel.cond(1) headmodel.cond(1) headmodel.cond(1) headmodel.cond(1)];
  case 2
    headmodel.r    = [headmodel.r(1)    headmodel.r(2)    headmodel.r(2)    headmodel.r(2)];
    headmodel.cond = [headmodel.cond(1) headmodel.cond(2) headmodel.cond(2) headmodel.cond(2)];
  case 3
    headmodel.r    = [headmodel.r(1)    headmodel.r(2)    headmodel.r(3)    headmodel.r(3)];
    headmodel.cond = [headmodel.cond(1) headmodel.cond(2) headmodel.cond(3) headmodel.cond(3)];
  case 4
    headmodel.r    = [headmodel.r(1)    headmodel.r(2)    headmodel.r(3)    headmodel.r(4)];
    headmodel.cond = [headmodel.cond(1) headmodel.cond(2) headmodel.cond(3) headmodel.cond(4)];
  otherwise
    error('not more than 4 spheres are supported');
end

% precompute the parameters for the series expansion
headmodel.t = eeg_leadfield4_prepare(headmodel, order);
