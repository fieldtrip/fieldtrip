function vol = ft_headmodel_infinite(varargin)

% FT_HEADMODEL_INFINITE returns an infinitely large homogenous
% volume conduction model. For EEG the volume conductor can be used
% to compute the leadfield of electric current dipoles, for MEG it
% can be used for computing the leadfield of magnmetic dipoles.
% 
% Use as
%   vol = ft_headmodel_infinite;
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

sourcemodel = ft_getopt(varargin, 'sourcemodel');

% this is an easy one
vol = [];

switch sourcemodel
case 'monopole'
  vol.type = 'infinite_monopole';  
case 'magneticdipole'
  vol.type = 'infinite_magneticdipole';  
case 'currentdipole'
  vol.type = 'infinite_currentdipole';  
otherwise
  vol.type = 'infinite';
end % switch

