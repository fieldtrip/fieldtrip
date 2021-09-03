function [pos_norm, chan_norm] = warp_mni(cfg, elec)

% WARP_MNI projects electrode channels to the MNI template brain.
% This volume-based registration technique can be used for the spatial
% normalization of electrodes based on a normalized anatomical volume. 
% To perform volume-based normalization, you first need to process the 
% subject's MRI with ft_volumenormalise.
%
% The configuration must contain the following options
%   cfg.volume         = anatomical volume, normalized with
%                        ft_volumenormalise
%
% See also FT_ELECTRODEREALIGN, FT_PREPARE_MESH

% Copyright (C) 2021, Arjen Stolk
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

% initial warp performed by ft_convert_coordsys
pos_norm  = ft_warp_apply(cfg.volume.initial, elec.elecpos, 'homogenous');
chan_norm = ft_warp_apply(cfg.volume.initial, elec.chanpos, 'homogenous');

% final warp performed by ft_volumenormalise
pos_norm  = ft_warp_apply(cfg.volume.params, pos_norm, 'individual2sn');
chan_norm = ft_warp_apply(cfg.volume.params, chan_norm, 'individual2sn');
