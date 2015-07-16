function headmodel = ft_headmodel_localspheres(geometry, grad, varargin)

% FT_HEADMODEL_LOCALSPHERES constructs a MEG volume conduction model in
% with a local sphere fitted to the head or brain surface for each separate
% channel
%
% This implements
%   Huang MX, Mosher JC, Leahy RM. "A sensor-weighted overlapping-sphere
%   head model and exhaustive head model comparison for MEG." Phys Med
%   Biol. 1999 Feb;44(2):423-40
%
% Use as
%   headmodel = ft_headmodel_localspheres(geom, grad, ...)
%
% Optional arguments should be specified in key-value pairs and can include
%   radius    = number, radius of sphere within which headshape points will
%               be included for the fitting algorithm
%   maxradius = number, if for a given sensor the fitted radius exceeds
%               this value, the radius and origin will be replaced with the
%               single sphere fit
%   baseline  = number
%   feedback  = boolean, true or false
%
% See also FT_PREPARE_HEADMODEL, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

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

% get the additional inputs and set the defaults
unit          = ft_getopt(varargin, 'unit');
feedback      = ft_getopt(varargin, 'feedback', true);
singlesphere  = ft_getopt(varargin, 'singlesphere', 'no');

if any(strcmp(varargin(1:2:end), 'unit')) || any(strcmp(varargin(1:2:end), 'units'))
  % the geometrical units should be specified in the input geometry
  error('the ''unit'' option is not supported any more');
end

% convert from 'yes'/'no' string into boolean value
feedback = istrue(feedback);

if isnumeric(geometry) && size(geometry,2)==3
  % assume that it is a Nx3 array with vertices
  % convert it to a structure, this is needed to determine the units further down
  geometry = struct('pnt', geometry);
elseif isstruct(geometry) && isfield(geometry,'bnd')
  % take the triangulated surfaces from the input structure
  geometry = geometry.bnd;
elseif ~(isstruct(geometry) && isfield(geometry,'pnt'))
  error('the input geometry should be a set of points or a single triangulated surface')
end

if isstruct(geometry) && numel(geometry)>1
  error('There must be only 1 geometry given as input');
end

% start with an empty volume conductor
headmodel = [];

% ensure that the geometry has units, estimate them if needed
geometry = ft_convert_units(geometry);

% ensure that it has consistent units
grad = ft_convert_units(grad, geometry.unit);

% copy the geometrical units into the volume conductor
headmodel.unit = geometry.unit;

% ensure that all defaults have the same user-defined units
radius    = ft_getopt(varargin, 'radius',    scalingfactor('cm', headmodel.unit) * 8.5);
maxradius = ft_getopt(varargin, 'maxradius', scalingfactor('cm', headmodel.unit) * 20);
baseline  = ft_getopt(varargin, 'baseline',  scalingfactor('cm', headmodel.unit) * 5);

% get the points from the triangulated surface
geometry = geometry.pnt;

Nshape = size(geometry,1);
Nchan  = numel(grad.label);

% set up an empty figure
if istrue(feedback)
  clf
  hold on
  axis equal
  axis vis3d
  axis off
  drawnow
end

% plot all channels and headshape points
if istrue(feedback)
  cla
  ft_plot_sens(grad);
  ft_plot_mesh(geometry, 'vertexcolor', 'g', 'facecolor', 'none', 'edgecolor', 'none');
  drawnow
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(geometry);
fprintf('single sphere,   %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', Nshape, single_o(1), single_o(2), single_o(3), single_r);

headmodel = [];
if strcmp(singlesphere, 'yes')
  % only return a single sphere
  headmodel.r = single_r;
  headmodel.o = single_o;
  return;
end

% allocate empty matrices that will hold the results
headmodel.r = zeros(Nchan,1);    % radius of every sphere
headmodel.o = zeros(Nchan,3);    % origin of every sphere
headmodel.label = cell(Nchan,1); % corresponding gradiometer channel label for every sphere

for chan=1:Nchan
  coilsel = find(grad.tra(chan,:)~=0);
  allpnt  = grad.coilpos(coilsel, :);   % position of all coils belonging to this channel
  allori  = grad.coilori(coilsel, :);   % orientation of all coils belonging to this channel
  
  if istrue(feedback)
    cla
    plot3(grad.coilpos(:,1), grad.coilpos(:,2), grad.coilpos(:,3), 'b.');   % all coils
    plot3(      allpnt(:,1),       allpnt(:,2),       allpnt(:,3), 'r*');     % this channel in red
  end
  
  % determine the average position and orientation of this channel
  thispnt = mean(allpnt,1);
  [u, s, v] = svd(allori);
  thisori = v(:,1)';
  if dot(thispnt,thisori)<0
    % the orientation should be outwards pointing
    thisori = -thisori;
  end
  
  % compute the distance from every coil along this channels orientation
  dist = zeros(size(coilsel));
  for i=1:length(coilsel)
    dist(i) = dot((allpnt(i,:)-thispnt), thisori);
  end
  
  [m, i] = min(dist);
  % check whether the minimum difference is larger than a typical distance
  if abs(m)>(baseline/4)
    % replace the position of this channel by the coil that is the closest to the head (axial gradiometer)
    % except when the center of the channel is approximately just as good (planar gradiometer)
    thispnt = allpnt(i,:);
  end
  
  % find the headshape points that are close to this channel
  dist = sqrt(sum((geometry-repmat(thispnt,Nshape,1)).^2, 2));
  shapesel = find(dist<radius);
  if feedback
    ft_plot_mesh(geometry(shapesel,:), 'vertexcolor', 'g');
    drawnow
  end
  
  % fit a sphere to these headshape points
  if length(shapesel)>10
    [o, r] = fitsphere(geometry(shapesel,:));
    fprintf('channel = %s, %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', grad.label{chan}, length(shapesel), o(1), o(2), o(3), r);
  else
    fprintf('channel = %s, not enough surface points, using all points\n', grad.label{chan});
    o = single_o;
    r = single_r;
  end
  
  if r > maxradius
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
