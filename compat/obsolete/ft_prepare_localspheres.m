function [headmodel, cfg] = ft_prepare_localspheres(cfg, mri)

% FT_PREPARE_LOCALSPHERES is deprecated, please use FT_PREPARE_HEADMODEL and
% FT_PREPARE_MESH
%
% See also FT_PREPARE_HEADMODEL

% Copyright (C) 2005-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
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

warning('FT_PREPARE_LOCALSPHERES is deprecated, please use FT_PREPARE_HEADMODEL with cfg.method = ''localspheres'' instead.')

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar mri
ft_preamble provenance mri


% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'spheremesh', 'numvertices'});

% set the defaults
if ~isfield(cfg, 'radius'),        cfg.radius = 8.5;        end
if ~isfield(cfg, 'maxradius'),     cfg.maxradius = 20;      end
if ~isfield(cfg, 'baseline'),      cfg.baseline = 5;        end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'yes';    end
if ~isfield(cfg, 'smooth');        cfg.smooth    = 5;       end % in voxels
if ~isfield(cfg, 'threshold'),     cfg.threshold = 0.5;     end % relative
if ~isfield(cfg, 'numvertices'),   cfg.numvertices = [];    end
if ~isfield(cfg, 'singlesphere'),  cfg.singlesphere = 'no'; end
if ~isfield(cfg, 'headshape'),     cfg.headshape = [];      end

hasmri = exist('mri', 'var'); % note that nargin will not work in case of cfg.inputfile
if hasmri
  headshape = ft_prepare_mesh(cfg, mri);
elseif isfield(cfg,'headshape') && nargin == 1
  if ischar(cfg.headshape)
    headshape = ft_read_headshape(cfg.headshape);
  else
    headshape = cfg.headshape;
  end
else
  ft_error('no head shape available')
end

% read the gradiometer definition from file or copy it from the configuration
if isfield(cfg, 'gradfile')
  grad = ft_read_sens(cfg.gradfile, 'senstype', 'meg');
else
  grad = cfg.grad;
end

Nshape = size(headshape.pos, 1);
Nchan  = size(grad.tra, 1);

% set up an empty figure
if strcmp(cfg.feedback, 'yes')
  clf
  hold on
  axis equal
  axis vis3d
  axis off
  drawnow
end

% plot all channels and headshape points
if strcmp(cfg.feedback, 'yes')
  cla
  ft_plot_sens(grad);
  ft_plot_mesh(headshape, 'vertexcolor', 'g', 'facecolor', 'none', 'edgecolor', 'none');
  drawnow
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(headshape.pos);
fprintf('single sphere,   %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', Nshape, single_o(1), single_o(2), single_o(3), single_r);

headmodel = [];

if strcmp(cfg.singlesphere, 'yes')
  % only return a single sphere
  headmodel.r = single_r;
  headmodel.o = single_o;
  return;
end

% start with an empty structure that will hold the results
headmodel.r = zeros(Nchan,1);    % radius of every sphere
headmodel.o = zeros(Nchan,3);    % origin of every sphere
headmodel.label = cell(Nchan,1); % corresponding gradiometer channel label for every sphere

for chan=1:Nchan
  coilsel = find(grad.tra(chan,:)~=0);
  allpos  = grad.coilpos(coilsel, :);   % position of all coils belonging to this channel
  allori  = grad.coilori(coilsel, :);   % orientation of all coils belonging to this channel

  if strcmp(cfg.feedback, 'yes')
    cla
    plot3(grad.coilpos(:,1), grad.coilpos(:,2), grad.coilpos(:,3), 'b.');   % all coils
    plot3(allpos(:,1), allpos(:,2), allpos(:,3), 'r*');     % this channel in red
  end

  % determine the average position and orientation of this channel
  thispos = mean(allpos,1);
  [u, s, v] = svd(allori);
  thisori = v(:,1)';
  if dot(thispos,thisori)<0
    % the orientation should be outwards pointing
    thisori = -thisori;
  end

  % compute the distance from every coil along this channels orientation
  dist = zeros(size(coilsel));
  for i=1:length(coilsel)
    dist(i) = dot((allpos(i,:)-thispos), thisori);
  end

  [m, i] = min(dist);
  % check whether the minimum difference is larger than a typical distance
  if abs(m)>(cfg.baseline/4)
    % replace the position of this channel by the coil that is the closest to the head (axial gradiometer)
    % except when the center of the channel is approximately just as good (planar gradiometer)
    thispos = allpos(i,:);
  end

  % find the headshape points that are close to this channel
  dist = sqrt(sum((headshape.pos-repmat(thispos,Nshape,1)).^2, 2));
  shapesel = find(dist<cfg.radius);
  if strcmp(cfg.feedback, 'yes')
    ft_plot_mesh(headshape.pos(shapesel,:), 'vertexcolor', 'g');
    drawnow
  end

  % fit a sphere to these headshape points
  if length(shapesel)>10
    [o, r] = fitsphere(headshape.pos(shapesel,:));
    fprintf('channel = %s, %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', grad.label{chan}, length(shapesel), o(1), o(2), o(3), r);
  else
    fprintf('channel = %s, not enough surface points, using all points\n', grad.label{chan});
    o = single_o;
    r = single_r;
  end

  if r > cfg.maxradius
    fprintf('channel = %s, not enough surface points, using all points\n', grad.label{chan});
    o = single_o;
    r = single_r;
  end

  % add this sphere to the volume conductor
  headmodel.o(chan,:)   = o;
  headmodel.r(chan)     = r;
  headmodel.label{chan} = grad.label{chan};
end % for all channels

headmodel.type = 'localspheres';

% ensure that the geometrical units are specified
headmodel = ft_determine_units(headmodel);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug

ft_postamble previous mri
ft_postamble provenance headmodel
ft_postamble history    headmodel
