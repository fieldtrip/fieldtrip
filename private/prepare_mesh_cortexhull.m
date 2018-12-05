function headshape = prepare_mesh_cortexhull(cfg)

% PREPARE_MESH_CORTEXHULL creates a mesh representing the cortex hull, i.e.
% the smoothed envelope around the pial surface created by FreeSurfer
%
% This function relies on the FreeSurfer and iso2mesh software packages
%
% Configuration options:
%   cfg.headshape    = a filename containing the pial surface computed by
%                      FreeSurfer recon-all ('/path/to/surf/lh.pial')
%   cfg.fshome       = FreeSurfer folder location
%                      (default: '/Applications/freesurfer')
%   cfg.resolution   = (optional, default: 1) resolution of the volume
%                      delimited by headshape being floodfilled by mris_fill
%   cfg.outer_surface_sphere = (optional, default: 15) diameter of the sphere
%                      used by make_outer_surface to close the sulci using
%                      morphological operations.
%   cfg.smooth_steps = (optional) number of (shrinking) smoothing
%                      iterations (default: 60)
%   cfg.laplace_steps = (optional) number of additional (non-shrinking)
%                      smoothing iterations (default: 0)
%
% See also FT_PREPARE_MESH

% Copyright (C) 2012-2018, Arjen Stolk, Gio Piantoni, Andrew Dykstra
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


% get the default options
surf                 = ft_getopt(cfg, 'headshape');
fshome               = ft_getopt(cfg, 'fshome', '/Applications/freesurfer');
resolution           = ft_getopt(cfg, 'resolution', 1);
outer_surface_sphere = ft_getopt(cfg, 'outer_surface_sphere', 15);
smooth_steps         = ft_getopt(cfg, 'smooth_steps', 60);
correcthull          = ft_getopt(cfg, 'correcthull', true);
laplace_steps        = ft_getopt(cfg, 'laplace_steps', 0);

% add the FreeSurfer environment
fprintf('adding the FreeSurfer environment\n')
addpath([fshome '/matlab']); % where make_outer_surface is located
setenv('FREESURFER_HOME', fshome);
PATH = getenv('PATH');
setenv('PATH', [PATH ':' fshome '/bin']); % where mris_fill is located

% temporary files
surf_filled   = [tempname() '_pial.filled.mgz'];
surf_outer    = [tempname() '_pial_outer'];
surf_smooth   = [tempname() '_pial_smooth'];

% fill the mesh
cmd = sprintf('mris_fill -c -r %d %s %s', resolution, surf, surf_filled);
system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);

% make outer surface
make_outer_surface(surf_filled, outer_surface_sphere, surf_outer)

% smooth using mris_smooth (this shrinks the mesh a bit)
cmd = sprintf('mris_smooth -nw -n %d %s %s', smooth_steps, surf_outer, surf_smooth);
system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);

% expand the mesh using mris_expand (to compensate for shrinkage if needed)
if correcthull
  % quantify hull shrinkage
  surf_nosmooth = [tempname() '_pial_nosmooth'];
  cmd = sprintf('mris_smooth -nw -n %d %s %s', 0, surf_outer, ...
    surf_nosmooth);
  system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]); % FIXME: maybe one iter?
  smooth = ft_read_headshape(surf_smooth);
  nosmooth = ft_read_headshape(surf_nosmooth);
  center = mean(nosmooth.pos,1); % FIXME: alternatively used boundary_mesh here to quantify shrinkage (problem is that orig mesh has more vertices)
  dist_s = sqrt( (smooth.pos(:,1)-center(:,1)).^2 + (smooth.pos(:,2)-center(:,2)).^2 + (smooth.pos(:,3)-center(:,3)).^2 );
  dist_ns = sqrt( (nosmooth.pos(:,1)-center(:,1)).^2 + (nosmooth.pos(:,2)-center(:,2)).^2 + (nosmooth.pos(:,3)-center(:,3)).^2 );
  idx = dist_ns>dist_s; % find nodes extending past the smooth hull
  % correct for shrinkage if needed
  if any(idx)
    expansion = zeros(size(smooth.pos,1),1);
    expansion(idx) = sqrt( (nosmooth.pos(idx,1)-smooth.pos(idx,1)).^2 + (nosmooth.pos(idx,2)-smooth.pos(idx,2)).^2 + (nosmooth.pos(idx,3)-smooth.pos(idx,3)).^2 );
    write_curv([fileparts(tempname()) filesep 'expansion'], expansion, size(smooth.tri,1));
    surf_smooth2  = [tempname() '_pial_smooth2'];
    cmd = sprintf('mris_expand -thickness -thickness_name %s %s %d %s', [fileparts(tempname()) filesep 'expansion'], surf_smooth, -1.0, surf_smooth2);
    system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);
    headshape = ft_read_headshape(surf_smooth2); % FIXME: this output is expanded at undesired locations
  else
    fprintf('no hull correction needed\n');
  end
end

% smooth using iso2mesh (non-shrinking)
if laplace_steps >= 1
  ft_hastoolbox('iso2mesh',1);
  fprintf('non-shrinking smoothing for %d iterations\n', laplace_steps)
  conn = meshconn(headshape.tri, size(headshape.pos,1)); % determine neighbors
  headshape.pos = smoothsurf(headshape.pos, [], conn, laplace_steps, 0, 'laplacianhc', .2);
end
