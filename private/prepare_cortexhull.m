function headshape = prepare_cortexhull(cfg)

% PREPARE_CORTEXHULL creates a mesh representing the cortex hull, i.e. the
% smoothed envelope around the pial surface created by freesurfer.
% PREPARE_CORTEXHULL relies on freesurfer's command line functions and
% 'make_outer_surface' in the freesurfer/matlab folder
%
% Configuration options:
%   cfg.method      = 'cortexhull'
%   cfg.headshape   = a filename containing the pial surface computed by
%                     freesurfer recon-all ('/path/to/surf/lh.pial')
%   cfg.resolution  = (optional, default: 1) resolution of the volume
%                     delimited by headshape being floodfilled by mris_fill
%   cfg.outer_surface_sphere = (optional, default: 15) diameter of the sphere
%                     used by make_outer_surface to close the sulci using
%                     morphological operations.
%   cfg.smooth_steps = (optional, default: 60) number of smoothing iterations
%                     performed by mris_smooth
%
% Error that 'mris_fill' was not found means that freesurfer is not installed.
% Error that 'make_outer_surface' is not in matlab, means that you need to add
% 'freesurfer/matlab' to your path.
%
% See also FT_PREPARE_MESH

% Copyright (C) 2012-2016, Gio Piantoni, Andrew Dykstra
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

disp('Please cite: Dykstra et al. 2012 Neuroimage PMID: 22155045')

% get the default options
resolution = ft_getopt(cfg, 'resolution', 1);
outer_surface_sphere = ft_getopt(cfg, 'outer_surface_sphere', 15);
smooth_steps = ft_getopt(cfg, 'smooth_steps', 60);

headshape = ft_getopt(cfg, 'headshape');

% temporary files
surf_filled = [tempname() '_pial.filled.mgz'];
surf_outer = [tempname() '_pial_outer'];
surf_smooth = [tempname() '_pial_smooth'];

cmd = sprintf('mris_fill -c -r %d %s %s', resolution, headshape, ...
  surf_filled);
system(cmd);

make_outer_surface(surf_filled, outer_surface_sphere, surf_outer)

cmd = sprintf('mris_smooth -nw -n %d %s %s', smooth_steps, surf_outer, ...
  surf_smooth);
system(cmd);

headshape = ft_read_headshape(surf_smooth);
