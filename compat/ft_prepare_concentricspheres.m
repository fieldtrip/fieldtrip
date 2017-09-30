function [headmodel, cfg] = ft_prepare_concentricspheres(cfg)

% FT_PREPARE_CONCENTRICSPHERES is deprecated, please use FT_PREPARE_HEADMODEL and
% FT_PREPARE_MESH
%
% See also FT_PREPARE_HEADMODEL

% Copyright (C) 2009, Vladimir Litvak & Robert Oostenveld
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

warning('FT_PREPARE_CONCENTRICSPHERES is deprecated, please use FT_PREPARE_HEADMODEL with cfg.method = ''concentricspheres'' instead.')

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

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
    ft_warning('using default conductivity values');
    cfg.conductivity = 1;
  elseif numel(headshape)==3
    ft_warning('using default conductivity values');
    cfg.conductivity = [1 1/80 1]*0.33;
  else
    % for a 2 or 4 sphere model the order of the compartments is potentially ambiguous, hence no default should be supplied
    ft_error('a conductivity value should be specified for each compartment');
  end
end

if strcmp(cfg.fitind, 'all')
  fitind = 1:numel(headshape);
else
  fitind = cfg.fitind;
end

% concatenate the vertices of all surfaces
pos = [];
for i = fitind
  pos = [pos ; headshape(i).pos];
end
% remove double vertices
pos = unique(pos, 'rows');

Npos = size(pos, 1);

% set up an empty figure
if strcmp(cfg.feedback, 'yes')
  clf
  hold on
  axis equal
  axis vis3d
  axis off
  drawnow
  colors = {'b', 'y', 'm', 'r'};
  [sphere_pos, sphere_tri] = icosahedron162;
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(pos);
fprintf('initial sphere: number of surface points = %d\n', Npos);
fprintf('initial sphere: center = [%.1f %.1f %.1f]\n', single_o(1), single_o(2), single_o(3));
fprintf('initial sphere: radius = %.1f\n', single_r);

% fit the radius of each concentric sphere to the corresponding surface points
headmodel = [];
headmodel.o = single_o;
for i = 1:numel(headshape)
  dist     = sqrt(sum(((headshape(end-i+1).pos - repmat(single_o, size(headshape(end-i+1).pos,1), 1)).^2), 2));
  headmodel.r(i) = mean(dist);

  if strcmp(cfg.feedback, 'yes')
    if ~isfield(headshape(end-i+1), 'tri')
      headshape(end-i+1).tri = [];
    end

    % plot the original surface
    bndtmp = [];
    bndtmp.pos = headshape(end-i+1).pos;
    bndtmp.tri = headshape(end-i+1).tri;
    ft_plot_mesh(bndtmp,'facecolor','none')

    % plot the sphere surface
    bndtmp = [];
    bndtmp.pos = sphere_pos*headmodel.r(i) + repmat(single_o, size(sphere_pos, 1), 1);
    bndtmp.tri = sphere_tri;
    ft_plot_mesh(bndtmp,'edgecolor',colors{mod(i, numel(colors)) + 1},'facecolor','none');
  end
end

if numel(cfg.conductivity)==numel(headshape)
  headmodel.cond = cfg.conductivity;
else
  ft_error('incorrect specification of cfg.conductivity');
end

headmodel.type = 'concentricspheres';

% ensure that the geometrical units are specified
headmodel = ft_determine_units(headmodel);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance headmodel
ft_postamble history    headmodel
