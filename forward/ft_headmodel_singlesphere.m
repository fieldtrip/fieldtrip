function headmodel = ft_headmodel_singlesphere(mesh, varargin)

% FT_HEADMODEL_SINGLESPHERE creates a volume conduction model of the
% head by fitting a spherical model to a set of points that describe
% the head surface.
%
% For MEG this implements Cuffin BN, Cohen D.  "Magnetic fields of
% a dipole in special volume conductor shapes" IEEE Trans Biomed Eng.
% 1977 Jul;24(4):372-81.
%
% Use as
%   headmodel = ft_headmodel_singlesphere(mesh, ...)
%
% Optional arguments should be specified in key-value pairs and can include
%   conductivity     = number, conductivity of the sphere
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% FIXME document both EEG and MEG case

% Copyright (C) 2012-2013, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

% get the optional arguments
conductivity = ft_getopt(varargin, 'conductivity', 1);

if any(strcmp(varargin(1:2:end), 'unit')) || any(strcmp(varargin(1:2:end), 'units'))
  % the geometrical units should be specified in the input mesh
  ft_error('the ''unit'' option is not supported any more');
end

if isnumeric(mesh) && size(mesh,2)==3
  % assume that it is a Nx3 array with vertices
  % convert it to a structure, this is needed to determine the units further down
  mesh = struct('pos', mesh);
elseif isstruct(mesh) && isfield(mesh,'bnd')
  % take the triangulated surfaces from the input structure
  mesh = mesh.bnd;
end

% replace pnt with pos
mesh = fixpos(mesh);

if ~isstruct(mesh) || numel(mesh)>1 || ~isfield(mesh, 'pos')
  ft_error('the input mesh should be a set of points or a single triangulated surface')
end

if numel(conductivity)~=1
  ft_error('the conductivity should be a single number')
end

if numel(mesh)~=1
  ft_error('fitting a single sphere requires a single mesh')
end

% start with an empty volume conductor
headmodel = [];

% ensure that the mesh has units, estimate them if needed
mesh = ft_determine_units(mesh);

% copy the geometrical units into the volume conductor
headmodel.unit = mesh.unit;

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(mesh.pos);

headmodel.r    = single_r;
headmodel.o    = single_o;
headmodel.cond = conductivity;
headmodel.type = 'singlesphere';

fprintf('single sphere: radius = %.1f, conductivity = %f\n', headmodel.r, headmodel.cond);
