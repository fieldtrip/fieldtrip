function [pos, lraxis, centre, radius] = hartmut_eyemodel(eyecfg, coordsys, unit)

% HARTMUT_EYEMODEL returns candidate positions for ocular sources in one eye and the
% left-right symmetry axis of the coordinate system. It is used by FT_DIPOLEFITTING for
% the HArtMuT extension, which fits the eyes as a mirror-symmetric dipole pair.
%
% Use as
%   [pos, lraxis, centre, radius] = hartmut_eyemodel(eyecfg, coordsys, unit)
% where the input is
%   eyecfg   = structure with fields radius, interocular and offset (all in mm), and an
%              optional field pos with the centre of one eye in head coordinates
%   coordsys = string describing the coordinate system, see FT_DETERMINE_COORDSYS
%   unit     = string with the geometrical units, e.g. 'mm' or 'cm'
% and the output is
%   pos      = Nx3 matrix with candidate positions in one eye, empty when they cannot be determined
%   lraxis   = the left-right symmetry axis, 'x', 'y' or 'z', empty when it cannot be determined
%   centre   = 1x3 vector with the centre of one eye in head coordinates, empty when not determined
%   radius   = scalar radius of the eye region in head coordinates
%
% The candidate positions are returned for a single eye; the partner in the other eye is
% its mirror image across the midsagittal plane. The default eye geometry is taken from the HArtMuT
% NYhead (MNI ICBM152) model, available from https://github.com/harmening/HArtMuT/tree/main/HArtMuTmodels;
% the ICBM152 head has no eye cavities and its eyes sit higher than in the Colin27 head, so the
% defaults stay clear of the skull when transformed to an individual.
% The defaults therefore only apply to an MNI-like coordinate system. When neither the eye centre
% nor the symmetry axis can be determined, the caller should fall back to fitting a single dipole.
%
% See also FT_DIPOLEFITTING, COORDSYS2LABEL, FT_INSIDE_HEADMODEL

% Copyright (C) 2026, Nils Harmening
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

pos    = zeros(0,3);
lraxis = '';
name   = {'x', 'y', 'z'};

% determine the centre and radius of the left eye in head coordinates
scale  = ft_scalingfactor('mm', unit);
radius = eyecfg.radius * scale;
if isfield(eyecfg, 'pos') && ~isempty(eyecfg.pos)
  % use the eye centre specified by the user
  lefteye = eyecfg.pos(1,:);
elseif ismember(lower(coordsys), {'mni', 'spm', 'tal', 'acpc'})
  % use the default eye centre, which is in MNI coordinates; an MNI-like coordinate system
  % is RAS, so x is left-right, y is anterior and z is superior
  lefteye = [-eyecfg.interocular/2, eyecfg.offset(1), eyecfg.offset(2)] * scale;
else
  % the default MNI eye centre does not apply and the user did not specify one
  lefteye = [];
end
centre = lefteye;

% determine the left-right axis, this is the axis whose two directions are left and right
[labelx, labely, labelz] = coordsys2label(coordsys, 1, true);
labels = {labelx, labely, labelz};
for i=1:3
  if all(ismember(labels{i}, {'left', 'right'}))
    lraxis = name{i};
  end
end

% both the eye centre and the symmetry axis are required to mirror the eye onto the other
% side, otherwise the caller should fall back to fitting a single dipole
if isempty(lefteye) || isempty(lraxis)
  return
end

% construct a regular grid of candidate positions inside the spherical eye region
step = radius/2;
[gx, gy, gz] = ndgrid(-radius:step:radius, -radius:step:radius, -radius:step:radius);
grid = [gx(:) gy(:) gz(:)];
grid = grid(sum(grid.^2, 2)<=radius^2, :); % keep the positions inside the sphere
pos  = grid + lefteye;
