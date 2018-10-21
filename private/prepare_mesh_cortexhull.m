function headshape = prepare_mesh_cortexhull(cfg)

% PREPARE_MESH_CORTEXHULL creates a mesh representing the cortex hull, i.e. the
% smoothed envelope around the pial surface created by FreeSurfer.
%
% This function relies on FreeSurfer's command line functions and
% 'make_outer_surface' in the FreeSurfer/matlab folder
%
% Configuration options:
%   cfg.method       = 'cortexhull'
%   cfg.headshape    = a filename containing the pial surface computed by
%                     FreeSurfer recon-all ('/path/to/surf/lh.pial')
%   cfg.resolution   = (optional, default: 1) resolution of the volume
%                     delimited by headshape being floodfilled by mris_fill
%   cfg.fshome       = FreeSurfer folder location (default: '/Applications/freesurfer')
%   cfg.outer_surface_sphere = (optional, default: 15) diameter of the sphere
%                     used by make_outer_surface to close the sulci using
%                     morphological operations.
%   cfg.smooth_steps = (optional, default: 60) number of smoothing iterations
%                     performed by mris_smooth
%
% See also FT_PREPARE_MESH

% Copyright (C) 2012-2016, Gio Piantoni, Andrew Dykstra, Arjen Stolk
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
resolution           = ft_getopt(cfg, 'resolution', 1);
outer_surface_sphere = ft_getopt(cfg, 'outer_surface_sphere', 15);
smooth_steps         = ft_getopt(cfg, 'smooth_steps', 60);
headshape            = ft_getopt(cfg, 'headshape');
fshome               = ft_getopt(cfg, 'fshome', '/Applications/freesurfer');

% add the FreeSurfer environment 
fprintf('adding the FreeSurfer environment\n')
addpath([fshome '/matlab']); % where make_outer_surface is located
setenv('FREESURFER_HOME', fshome);
PATH = getenv('PATH');
setenv('PATH', [PATH ':' fshome '/bin']); % where mris_fill is located

% temporary files
surf_filled  = [tempname() '_pial.filled.mgz'];
surf_outer   = [tempname() '_pial_outer'];
surf_smooth  = [tempname() '_pial_smooth'];

cmd = sprintf('mris_fill -c -r %d %s %s', resolution, headshape, ...
  surf_filled);
system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);

make_outer_surface(surf_filled, outer_surface_sphere, surf_outer)

cmd = sprintf('mris_smooth -nw -n %d %s %s', smooth_steps, surf_outer, ...
  surf_smooth);
system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);

headshape = ft_read_headshape(surf_smooth);
