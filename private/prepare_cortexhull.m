function headshape = prepare_cortexhull(cfg)

% FT_PREPARE_MESH creates a triangulated surface mesh for the volume
% conduction model. The mesh can either be selected manually from raw
% mri data or can be generated starting from a segmented volume
% information stored in the mri structure. The result is a bnd
% structure which contains the information about all segmented surfaces
% related to mri and are expressed in world coordinates.
%
% Use as
%   bnd = ft_prepare_mesh(cfg, mri)
%   bnd = ft_prepare_mesh(cfg, seg)
%
% Configuration options:
%   cfg.method      = string, can be 'interactive', 'projectmesh', 'iso2mesh', 'isosurface',
%                     'headshape', 'hexahedral', 'tetrahedral', 'cortexhull'
%   cfg.tissue      = cell-array with tissue types or numeric vector with integer values
%   cfg.numvertices = numeric vector, should have same number of elements as cfg.tissue
%   cfg.downsample  = integer number (default = 1, i.e. no downsampling), see FT_VOLUMEDOWNSAMPLE
%   cfg.headshape   = (optional) a filename containing headshape, a Nx3 matrix with surface
%                     points, or a structure with a single or multiple boundaries


% obligatory fields
% cfg.headshape: path to lh.pial or rh.pial

% Original author: Andrew Dykstra
% adapted: Gio Piantoni (gio@gpiantoni.com)

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
