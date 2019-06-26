function headmodel = ft_headmodel_localspheres(mesh, grad, varargin)

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
%   headmodel = ft_headmodel_localspheres(mesh, grad, ...)
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

% get the additional inputs and set the defaults
feedback      = ft_getopt(varargin, 'feedback', true);
singlesphere  = ft_getopt(varargin, 'singlesphere', 'no');

if any(strcmp(varargin(1:2:end), 'unit')) || any(strcmp(varargin(1:2:end), 'units'))
  % the geometrical units should be specified in the input mesh
  ft_error('the ''unit'' option is not supported any more');
end

% convert from 'yes'/'no' string into boolean value
feedback = istrue(feedback);

if isnumeric(mesh) && size(mesh,2)==3
  % assume that it is a Nx3 array with vertices
  % convert it to a structure, this is needed to determine the units further down
  mesh = struct('pos', mesh);
elseif isstruct(mesh) && isfield(mesh, 'bnd')
  % take the triangulated surfaces from the input structure
  mesh = mesh.bnd;
end

% replace pnt with pos
mesh = fixpos(mesh);

if ~isstruct(mesh) || numel(mesh)>1 || ~isfield(mesh, 'pos')
  ft_error('the input mesh should be a set of points or a single triangulated surface')
end

if isstruct(mesh) && numel(mesh)>1
  ft_error('There must be only 1 mesh given as input');
end

% start with an empty volume conductor
headmodel = [];

% ensure that the mesh has units, estimate them if needed
mesh = ft_determine_units(mesh);

% ensure that it has a consistent representation and consistent units
grad = ft_datatype_sens(grad);
grad = ft_convert_units(grad, mesh.unit);

% copy the geometrical units into the volume conductor
headmodel.unit = mesh.unit;

% ensure that all defaults have the same user-defined units
radius    = ft_getopt(varargin, 'radius',    ft_scalingfactor('cm', headmodel.unit) * 8.5);
maxradius = ft_getopt(varargin, 'maxradius', ft_scalingfactor('cm', headmodel.unit) * 20);
baseline  = ft_getopt(varargin, 'baseline',  ft_scalingfactor('cm', headmodel.unit) * 5);

% get the points from the triangulated surface
mesh = mesh.pos;

Nshape = size(mesh,1);
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
  ft_plot_mesh(mesh, 'vertexcolor', 'g', 'facecolor', 'none', 'edgecolor', 'none');
  drawnow
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(mesh);
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

allpos = grad.chanpos;
allori = grad.chanori;
if any(~isfinite(allpos(:)))
  % probably the data has been backprojected from a component structure,
  % use chanposold, if possible
  if isfield(grad, 'chanposold')
    allpos = grad.chanposold;
    allori = grad.chanoriold;
  else
    ft_error('cannot extract channel position information from sensor description');
  end
end

for chan=1:Nchan
  thispos  = allpos(chan,:);
  thisori  = allori(chan,:);
 
  if istrue(feedback)
    cla
    plot3( allpos(:,1),  allpos(:,2),  allpos(:,3), 'b.');   % all channels
    plot3(thispos(:,1), thispos(:,2), thispos(:,3), 'r*');     % this channel in red
  end
  
  if dot(thispos,thisori)<0
    % the orientation should be outwards pointing
    thisori = -thisori;
  end
    
  % find the headshape points that are close to this channel
  dist = sqrt(sum((mesh-repmat(thispos,Nshape,1)).^2, 2));
  shapesel = find(dist<radius);
  if feedback
    ft_plot_mesh(mesh(shapesel,:), 'vertexcolor', 'g');
    drawnow
  end
  
  % fit a sphere to these headshape points
  if length(shapesel)>10
    [o, r] = fitsphere(mesh(shapesel,:));
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
