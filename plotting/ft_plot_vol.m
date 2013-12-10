function ft_plot_vol(vol, varargin)

% FT_PLOT_VOL visualizes the boundaries in the vol structure constituting the
% geometrical information of the forward model
%
% Use as
%   hs = ft_plot_vol(vol, varargin)
%
% Graphic facilities are available for vertices, edges and faces. A list of
% the arguments is given below with the correspondent admitted choices.
%
%     'facecolor'     [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'vertexcolor'   [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'edgecolor'     [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'facealpha'     number between 0 and 1
%     'faceindex'     true or false
%     'vertexindex'   true or false
%
% Example
%   vol.r = [86 88 92 100];
%   vol.o = [0 0 40];
%   figure, ft_plot_vol(vol)

% Copyright (C) 2009, Cristiano Micheli
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

ws = warning('on', 'MATLAB:divideByZero');

% ensure that the volume conduction model description is up-to-date (Dec 2012)
vol = ft_datatype_headmodel(vol);

% get the optional input arguments
faceindex   = ft_getopt(varargin, 'faceindex', 'none');
vertexindex = ft_getopt(varargin, 'vertexindex', 'none');
vertexsize  = ft_getopt(varargin, 'vertexsize', 10);
facecolor   = ft_getopt(varargin, 'facecolor', 'white');
vertexcolor = ft_getopt(varargin, 'vertexcolor', 'none');
edgecolor   = ft_getopt(varargin, 'edgecolor'); % the default for this is set below
facealpha   = ft_getopt(varargin, 'facealpha', 1);
surfaceonly = ft_getopt(varargin, 'surfaceonly');

faceindex   = istrue(faceindex);   % yes=view the face number
vertexindex = istrue(vertexindex); % yes=view the vertex number

% we will probably need a sphere, so let's prepare one
[pnt, tri] = icosahedron2562;

% prepare a single or multiple triangulated boundaries
switch ft_voltype(vol)
  case {'singlesphere' 'concentricspheres'}
    vol.r = sort(vol.r);
    bnd = [];
    for i=1:length(vol.r)
      bnd(i).pnt(:,1) = pnt(:,1)*vol.r(i) + vol.o(1);
      bnd(i).pnt(:,2) = pnt(:,2)*vol.r(i) + vol.o(2);
      bnd(i).pnt(:,3) = pnt(:,3)*vol.r(i) + vol.o(3);
      bnd(i).tri = tri;
    end
    if isempty(edgecolor)
      edgecolor = 'none';
    end
    
  case 'localspheres'
    bnd = [];
    for i=1:length(vol.label)
      bnd(i).pnt(:,1) = pnt(:,1)*vol.r(i) + vol.o(i,1);
      bnd(i).pnt(:,2) = pnt(:,2)*vol.r(i) + vol.o(i,2);
      bnd(i).pnt(:,3) = pnt(:,3)*vol.r(i) + vol.o(i,3);
      bnd(i).tri = tri;
    end
    if isempty(edgecolor)
      edgecolor = 'none';
    end
    
  case {'bem', 'dipoli', 'asa', 'bemcp', 'singleshell' 'openmeeg'}
    % these already contain one or multiple triangulated surfaces for the boundaries
    bnd = vol.bnd;
    
  case 'simbio'
    % the ft_plot_mesh function below wants the SIMBIO tetrahedral or hexahedral mesh
    bnd = vol;
    
    % only plot the outer surface of the volume
    surfaceonly = true;
    
  case 'interpolate'
    xgrid = 1:vol.dim(1);
    ygrid = 1:vol.dim(2);
    zgrid = 1:vol.dim(3);
    [x, y, z] = ndgrid(xgrid, ygrid, zgrid);
    gridpos = ft_warp_apply(vol.transform, [x(:) y(:) z(:)]);
    
    % plot the dipole positions that are inside
    plot3(gridpos(vol.inside, 1), gridpos(vol.inside, 2), gridpos(vol.inside, 3), 'k.');
    
    % there is no boundary to be displayed
    bnd = [];
    
  case {'infinite' 'infinite_monopole' 'infinite_currentdipole' 'infinite_magneticdipole'}
    warning('there is nothing to plot for an infinite volume conductor')
    
    % there is no boundary to be displayed
    bnd = [];
    
  otherwise
    error('unsupported voltype')
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
warning(ws);
