function [vol, cfg] = ft_prepare_singleshell(cfg, mri)

% FT_PREPARE_SINGLESHELL creates a simple and fast method for the MEG forward
% calculation for one shell of arbitrary shape. This is based on a
% correction of the lead field for a spherical volume conductor by a
% superposition of basis functions, gradients of harmonic functions
% constructed from spherical harmonics.
%
% Use as
%   vol = ft_prepare_singleshell(cfg, seg), or
%   vol = ft_prepare_singleshell(cfg, mri), or
%   vol = ft_prepare_singleshell(cfg)
%
% If you do not use a segmented MRI, the configuration should contain
%   cfg.headshape   = a filename containing headshape, a structure containing a
%                     single triangulated boundary, or a Nx3 matrix with surface
%                     points
%   cfg.numvertices = number, to retriangulate the mesh with a sphere (default = 3000)
%                     instead of specifying a number, you can specify 'same' to keep the
%                     vertices of the mesh identical to the original headshape points
%
% The following options are relevant if you use a segmented MRI
%   cfg.smooth      = 'no' or the FWHM of the gaussian kernel in voxels (default = 5)
%   cfg.sourceunits = 'mm' or 'cm' (default is 'cm')
%   cfg.threshold   = 0.5, relative to the maximum value in the segmentation
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'mri',
% corresponding to the input structure.
%
% This function implements
%   G. Nolte, "The magnetic lead field theorem in the quasi-static
%   approximation and its use for magnetoencephalography forward calculation
%   in realistic volume conductors", Phys Med Biol. 2003 Nov 21;48(22):3637-52.
%
% See also FT_PREPARE_CONCENTRICSPHERES, FT_PREPARE_LOCALSPHERES,
% FT_PREPARE_BEMMODEL, FT_PREPARE_LEADFIELD, FT_PREPARE_MESH,
% FT_PREPARE_MESH_NEW

% TODO the spheremesh option should be renamed consistently with other mesh generation cfgs
% TODO shape should contain pnt as subfield and not be equal to pnt (for consistency with other use of shape)

% Copyright (C) 2006-2007, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar mri

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'spheremesh', 'numvertices'});
cfg = ft_checkconfig(cfg, 'deprecated', 'mriunits');

% set the defaults
if ~isfield(cfg, 'smooth');        cfg.smooth = 5;          end % in voxels
if ~isfield(cfg, 'sourceunits'),   cfg.sourceunits = 'cm';  end
if ~isfield(cfg, 'threshold'),     cfg.threshold = 0.5;     end % relative
if ~isfield(cfg, 'numvertices'),   cfg.numvertices = [];    end % approximate number of vertices in sphere

% the data is specified as input variable or input file
hasmri = exist('mri', 'var');

if hasmri && isempty(cfg.numvertices)
  cfg.numvertices = 3000;
end

if hasmri
  vol.bnd = ft_prepare_mesh(cfg, mri);
else
  vol.bnd = ft_prepare_mesh(cfg);
end

vol.type = 'singleshell';

% ensure that the geometrical units are specified
vol = ft_convert_units(vol);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
if hasmri
  ft_postamble previous mri
end
ft_postamble history vol

