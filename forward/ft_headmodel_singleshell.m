function headmodel = ft_headmodel_singleshell(mesh, varargin)

% FT_HEADMODEL_SINGLESHELL creates a volume conduction model of the
% head for MEG based on a realistic shaped surface of the inside of
% the skull.
% 
% The method implemented in this function allows for a simple and
% fast method for the MEG forward calculation for one shell of arbitrary
% shape, based on a correction of the lead field for a spherical
% volume conductor by a superposition of basis functions, gradients
% of harmonic functions constructed from spherical harmonics.
% 
% This function implements
%   G. Nolte, "The magnetic lead field theorem in the quasi-static
%   approximation and its use for magnetoencephalography forward calculation
%   in realistic volume conductors", Phys Med Biol. 2003 Nov 21;48(22):3637-52.
% 
% Use as
%   headmodel = ft_headmodel_singleshell(mesh, ...)
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

if ~isstruct(mesh) || ~isfield(mesh, 'pos')
  ft_error('the input mesh should be a set of points or a single triangulated surface')
end

% represent the mesh in a headmodel strucure
% the computational parameters will be added later on by ft_prepare_vol_sens
headmodel      = [];
headmodel.bnd  = mesh;
headmodel.type = 'singleshell';
headmodel = ft_determine_units(headmodel);

