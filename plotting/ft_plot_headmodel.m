function ft_plot_headmodel(headmodel, varargin)

% FT_PLOT_HEADMODEL visualizes the boundaries in the volume conduction model of the head as
% specified in the headmodel structure
%
% Use as
%   hs = ft_plot_headmodel(headmodel, varargin)
%
% Optional arguments should come in key-value pairs and can include
%   'facecolor'    = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r', or an Nx3 or Nx1 array where N is the number of faces
%   'vertexcolor'  = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r', or an Nx3 or Nx1 array where N is the number of vertices
%   'edgecolor'    = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%   'faceindex'    = true or false
%   'vertexindex'  = true or false
%   'facealpha'    = transparency, between 0 and 1 (default = 1)
%   'edgealpha'    = transparency, between 0 and 1 (default = 1)
%   'surfaceonly'  = true or false, plot only the outer surface of a hexahedral or tetrahedral mesh (default = false)
%   'unit'         = string, convert to the specified geometrical units (default = [])
%   'grad'         = gradiometer array, used in combination with local spheres model
%
% Example
%   headmodel   = [];
%   headmodel.r = [86 88 92 100];
%   headmodel.o = [0 0 40];
%   figure, ft_plot_headmodel(headmodel)
%
% See also FT_PREPARE_HEADMODEL, FT_PLOT_MESH, FT_PLOT_SENS

% Copyright (C) 2009, Cristiano Micheli
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

ws = ft_warning('on', 'MATLAB:divideByZero');

% ensure that the volume conduction model description is up-to-date (Dec 2012)
headmodel = ft_datatype_headmodel(headmodel);

% get the optional input arguments
faceindex   = ft_getopt(varargin, 'faceindex', 'none');
vertexindex = ft_getopt(varargin, 'vertexindex', 'none');
vertexsize  = ft_getopt(varargin, 'vertexsize', 10);
facecolor   = ft_getopt(varargin, 'facecolor', 'white');
vertexcolor = ft_getopt(varargin, 'vertexcolor', 'none');
edgecolor   = ft_getopt(varargin, 'edgecolor'); % the default for this is set below
facealpha   = ft_getopt(varargin, 'facealpha', 1);
surfaceonly = ft_getopt(varargin, 'surfaceonly');
unit        = ft_getopt(varargin, 'unit');
grad        = ft_getopt(varargin, 'grad');

if ~isempty(unit)
  headmodel = ft_convert_units(headmodel, unit);
  if ~isempty(grad)
    grad = ft_convert_units(grad, unit);
  end
end

faceindex   = istrue(faceindex);   % yes=view the face number
vertexindex = istrue(vertexindex); % yes=view the vertex number

% we will probably need a sphere, so let's prepare one
[pos, tri] = mesh_sphere(2562);

% prepare a single or multiple triangulated boundaries
switch ft_headmodeltype(headmodel)
  case {'singlesphere' 'concentricspheres'}
    headmodel.r = sort(headmodel.r);
    bnd = repmat(struct(), numel(headmodel.r));
    for i=1:numel(headmodel.r)
      bnd(i).pos(:,1) = pos(:,1)*headmodel.r(i) + headmodel.o(1);
      bnd(i).pos(:,2) = pos(:,2)*headmodel.r(i) + headmodel.o(2);
      bnd(i).pos(:,3) = pos(:,3)*headmodel.r(i) + headmodel.o(3);
      bnd(i).tri = tri;
    end
    if isempty(edgecolor)
      edgecolor = 'none';
    end

  case 'localspheres'
    if ~isempty(grad)
      ft_notice('estimating point on head surface for each gradiometer');
      [headmodel, grad] = ft_prepare_vol_sens(headmodel, grad);
      [bnd.pos, bnd.tri] = headsurface(headmodel, grad);
    else
      ft_notice('plotting sphere for each gradiometer');
      bnd = repmat(struct(), numel(headmodel.label));
      for i=1:numel(headmodel.label)
        bnd(i).pos(:,1) = pos(:,1)*headmodel.r(i) + headmodel.o(i,1);
        bnd(i).pos(:,2) = pos(:,2)*headmodel.r(i) + headmodel.o(i,2);
        bnd(i).pos(:,3) = pos(:,3)*headmodel.r(i) + headmodel.o(i,3);
        bnd(i).tri = tri;
      end
    end
    if isempty(edgecolor)
      edgecolor = 'none';
    end

  case {'bem', 'dipoli', 'asa', 'bemcp', 'singleshell' 'openmeeg'}
    % these already contain one or multiple triangulated surfaces for the boundaries
    bnd = headmodel.bnd;

  case 'simbio'
    % the ft_plot_mesh function below wants the SIMBIO tetrahedral or hexahedral mesh
    bnd = headmodel;

    % only plot the outer surface of the volume
    surfaceonly = true;

  case 'interpolate'
    xgrid = 1:headmodel.dim(1);
    ygrid = 1:headmodel.dim(2);
    zgrid = 1:headmodel.dim(3);
    [x, y, z] = ndgrid(xgrid, ygrid, zgrid);
    gridpos = ft_warp_apply(headmodel.transform, [x(:) y(:) z(:)]);

    % plot the dipole positions that are inside
    plot3(gridpos(headmodel.inside, 1), gridpos(headmodel.inside, 2), gridpos(headmodel.inside, 3), 'k.');

    % there is no boundary to be displayed
    bnd = [];

  case {'infinite' 'infinite_monopole' 'infinite_currentdipole' 'infinite_magneticdipole'}
    ft_warning('there is nothing to plot for an infinite volume conductor')

    % there is no boundary to be displayed
    bnd = [];

  otherwise
    ft_error('unsupported headmodel type')
end

% all models except for the spherical ones
if isempty(edgecolor)
  edgecolor = 'k';
end

% plot the triangulated surfaces of the volume conduction model
for i=1:length(bnd)
  ft_plot_mesh(bnd(i),'faceindex',faceindex,'vertexindex',vertexindex, ...
    'vertexsize',vertexsize,'facecolor',facecolor,'edgecolor',edgecolor, ...
    'vertexcolor',vertexcolor,'facealpha',facealpha, 'surfaceonly', surfaceonly);
end

% revert to original state
ft_warning(ws);
