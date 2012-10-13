function [vol, cfg] = ft_prepare_concentricspheres(cfg)

% FT_PREPARE_CONCENTRICSPHERES creates a EEG volume conductor model with
% multiple concentric spheres.
%
% Use as
%   vol = ft_prepare_concentricspheres(cfg)
%
% The input configuration should contain
%   cfg.headshape     = a filename containing headshape, a Nx3 matrix with surface 
%                       points, or a structure with a single or multiple boundaries
%   cfg.conductivity  = conductivity values for the model (default = [0.3300 1 0.0042 0.3300])
%   cfg.fitind        = indices of shapes to use for fitting the center (default = 'all')
%   cfg.nonlinear     = 'yes' or 'no' (default = 'yes')
%   cfg.feedback      = 'yes' or 'no' (default = 'yes')
%
% Example:
%
%   % first create 4 surfaces that represent the inner_skull_surface, csf, outer_skull_surface and skin_surface
%   radius = [86 88 92 100];
%   headshape = [];
%   for i=1:4
%     pnt = randn(100,3);
%     for j=1:size(pnt,1)
%       pnt(j,:) = pnt(j,:) ./ norm(pnt(j,:));
%     end
%     headshape(i).pnt = radius(i) .* pnt + 0.1*randn(size(pnt));
%   end
%
%   % then construct a volume conduction model of the head by fitting 4 concentric spheres
%   cfg = [];
%   cfg.headshape    = headshape;
%   cfg.conductivity = [0.3300 1 0.0042 0.3300]
%   [vol, cfg] = prepare_concentricspheres(cfg)
%
% See also FT_PREPARE_MESH, FT_PREPARE_HEADMODEL, FT_PREPARE_LEADFIELD

% Copyright (C) 2009, Vladimir Litvak & Robert Oostenveld
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
ft_preamble trackconfig
ft_preamble provenance

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', 'nonlinear');

% set the defaults
if ~isfield(cfg, 'fitind'),        cfg.fitind = 'all';             end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'yes';           end
if ~isfield(cfg, 'conductivity'),  cfg.conductivity = [];          end % this should be specified by the user
if ~isfield(cfg, 'numvertices'),   cfg.numvertices = 'same';       end

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested config-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

% get the surface describing the head shape
headshape = prepare_mesh_headshape(cfg);

if isempty(cfg.conductivity)
  if numel(headshape)==1
    warning('using default conductivity values');
    cfg.conductivity = 1;
  elseif numel(headshape)==3
    warning('using default conductivity values');
    cfg.conductivity = [1 1/80 1]*0.33;
  else
    % for a 2 or 4 sphere model the order of the compartments is potentially ambiguous, hence no default should be supplied
    error('a conductivity value should be specified for each compartment');
  end
end

if strcmp(cfg.fitind, 'all')
  fitind = 1:numel(headshape);
else
  fitind = cfg.fitind;
end

% concatenate the vertices of all surfaces
pnt = [];
for i = fitind
  pnt = [pnt ; headshape(i).pnt];
end
% remove double vertices
pnt = unique(pnt, 'rows');

Npnt = size(pnt, 1);

% set up an empty figure
if strcmp(cfg.feedback, 'yes')
  clf
  hold on
  axis equal
  axis vis3d
  axis off
  drawnow
  colors = {'b', 'y', 'm', 'r'};
  [sphere_pnt, sphere_tri] = icosahedron162;
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(pnt);
fprintf('initial sphere: number of surface points = %d\n', Npnt);
fprintf('initial sphere: center = [%.1f %.1f %.1f]\n', single_o(1), single_o(2), single_o(3));
fprintf('initial sphere: radius = %.1f\n', single_r);

% fit the radius of each concentric sphere to the corresponding surface points
vol = [];
vol.o = single_o;
for i = 1:numel(headshape)
  dist     = sqrt(sum(((headshape(end-i+1).pnt - repmat(single_o, size(headshape(end-i+1).pnt,1), 1)).^2), 2));
  vol.r(i) = mean(dist);

  if strcmp(cfg.feedback, 'yes')
    if ~isfield(headshape(end-i+1), 'tri')
      headshape(end-i+1).tri = [];
    end

    % plot the original surface
    bndtmp = [];
    bndtmp.pnt = headshape(end-i+1).pnt;
    bndtmp.tri = headshape(end-i+1).tri;
    ft_plot_mesh(bndtmp,'facecolor','none')

    % plot the sphere surface
    bndtmp = [];
    bndtmp.pnt = sphere_pnt*vol.r(i) + repmat(single_o, size(sphere_pnt, 1), 1);
    bndtmp.tri = sphere_tri;
    ft_plot_mesh(bndtmp,'edgecolor',colors{mod(i, numel(colors)) + 1},'facecolor','none');
  end
end

if numel(cfg.conductivity)==numel(headshape)
  vol.c = cfg.conductivity;
else
  error('incorrect specification of cfg.conductivity');
end

vol.type = 'concentricspheres';

% ensure that the geometrical units are specified
vol = ft_convert_units(vol);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble provenance
ft_postamble history vol

