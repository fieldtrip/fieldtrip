function vol = ft_headmodel_singleshell(geometry, varargin)

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
%   vol = ft_headmodel_singleshell(geom, ...)
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

% if it contains more than 1 shell it retunrs an error
if isfield(geometry,'pnt')
  if numel(geometry)>1
    error('no more than 1 shell at a time is allowed')
  end
else
  error('the input should be a boundary')
end

% represent the geometry in a headmodel strucure
% the computational parameters will be added later on by ft_prepare_vol_sens
vol      = [];
vol.bnd  = geometry;
vol.type = 'singleshell';
if ~isfield(vol, 'unit')
  vol = ft_convert_units(vol);
end

