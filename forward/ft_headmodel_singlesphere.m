function vol = ft_headmodel_singlesphere(geometry, varargin)

% FT_HEADMODEL_SINGLESPHERE creates a volume conduction model of the
% head by fitting a spherical model to a set of points that describe
% the head surface.
%
% For MEG this implements Cuffin BN, Cohen D.  "Magnetic fields of
% a dipole in special volume conductor shapes" IEEE Trans Biomed Eng.
% 1977 Jul;24(4):372-81.
%
% Use as
%   vol = ft_headmodel_singlesphere(pnt, ...)
%
% Optional arguments should be specified in key-value pairs and can include
%   conductivity     = number, conductivity of the sphere
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% FIXME document the EEG case

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

% get the optional arguments
conductivity = ft_getopt(varargin, 'conductivity', 1);
unit         = ft_getopt(varargin, 'unit');

if length(conductivity)~=1
  error('the conductivity should be a single number')
end

% start with an empty volume conductor
vol = [];

if ~isempty(unit)
  vol.unit = unit;                       % use the user-specified units for the output
else
  geometry = ft_convert_units(geometry); % ensure that it has units, estimate them if needed
  vol.unit = geometry.unit;              % copy the geometrical units into the volume conductor
end

if isnumeric(geometry) && size(geometry,2)==3
  % assume that it is a Nx3 array with vertices
elseif isstruct(geometry) && isfield(geometry,'pnt') && numel(geometry)==1
  % get the points from the triangulated surface
  geometry = geometry.pnt;
else
  error('the input geometry should be a set of points or a single triangulated surface')
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(geometry);

vol.r = single_r;
vol.o = single_o;
vol.c = conductivity;
vol.type = 'singlesphere';
vol      = ft_convert_units(vol); % ensure the object to have a unit

