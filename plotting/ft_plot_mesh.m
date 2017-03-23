function [hs] = ft_plot_mesh(mesh, varargin)

% FT_PLOT_MESH visualizes the information of a mesh contained in the first
% argument mesh. The boundary argument (mesh) typically contains two fields
% called .pos and .tri referring to the vertices and the triangulation of
% the mesh.
%
% Use as
%   ft_plot_mesh(mesh, ...)
% or if you only want to plot the 3-D vertices
%   ft_plot_mesh(pos, ...)
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
%   'vertexmarker' = character, e.g. '.', 'o' or 'x' (default = '.')
%   'vertexsize'   = scalar or vector with the size for each vertex (default = 10)
%   'unit'         = string, convert to the specified geometrical units (default = [])
%
% If you don't want the faces, edges or vertices to be plotted, you should specify the color as 'none'.
%
% Example
%   [pos, tri] = icosahedron162;
%   mesh.pos = pos;
%   mesh.tri = tri;
%   ft_plot_mesh(mesh, 'facecolor', 'skin', 'edgecolor', 'none')
%   camlight
%
% See also TRIMESH, PATCH

% Copyright (C) 2009-2015, Robert Oostenveld
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

ws = warning('on', 'MATLAB:divideByZero');

% rename pnt into pos
mesh = fixpos(mesh);

if ~isstruct(mesh) && isnumeric(mesh) && size(mesh,2)==3
  % the input seems like a list of points, convert into something that resembles a mesh
  mesh = struct('pos', mesh);
end

% the input is a structure, but might also be a struct-array
if numel(mesh)>1
  % plot each of the boundaries
  for i=1:numel(mesh)
    ft_plot_mesh(mesh(i), varargin{:})
  end
  return
end

% get the optional input arguments
vertexcolor  = ft_getopt(varargin, 'vertexcolor');
if isfield(mesh, 'tri') && size(mesh.tri,1)>10000
  facecolor    = ft_getopt(varargin, 'facecolor',   'cortex_light');
  edgecolor    = ft_getopt(varargin, 'edgecolor',   'none');
else
  facecolor    = ft_getopt(varargin, 'facecolor',   'white');
  edgecolor    = ft_getopt(varargin, 'edgecolor',   'k');
end
faceindex    = ft_getopt(varargin, 'faceindex',   false);
vertexindex  = ft_getopt(varargin, 'vertexindex', false);
vertexsize   = ft_getopt(varargin, 'vertexsize',  10);
vertexmarker = ft_getopt(varargin, 'vertexmarker', '.');
facealpha    = ft_getopt(varargin, 'facealpha',   1);
edgealpha    = ft_getopt(varargin, 'edgealpha',   1);
tag          = ft_getopt(varargin, 'tag',         '');
surfaceonly  = ft_getopt(varargin, 'surfaceonly');  % default is handled below
unit         = ft_getopt(varargin, 'unit');

haspos   = isfield(mesh, 'pos');  % vertices
hastri   = isfield(mesh, 'tri');  % triangles   as a Mx3 matrix with vertex indices
hastet   = isfield(mesh, 'tet');  % tetraheders as a Mx4 matrix with vertex indices
hashex   = isfield(mesh, 'hex');  % hexaheders  as a Mx8 matrix with vertex indices
hasline  = isfield(mesh, 'line'); % line segments in 3-D
haspoly  = isfield(mesh, 'poly'); % polynomial surfaces in 3-D
hascolor = isfield(mesh, 'color'); % color code for vertices

if hastet && isempty(surfaceonly)
  warning('only visualizing the outer surface of the tetrahedral mesh, see the "surfaceonly" option')
  surfaceonly = true;
elseif hashex && isempty(surfaceonly)
  warning('only visualizing the outer surface of the hexahedral mesh, see the "surfaceonly" option')
  surfaceonly = true;
else
  surfaceonly = false;
end

if ~isempty(unit)
  mesh = ft_convert_units(mesh, unit);
end

if surfaceonly
  mesh = mesh2edge(mesh);
  % update the flags that indicate which surface/volume elements are present
  hastri   = isfield(mesh, 'tri');  % triangles   as a Mx3 matrix with vertex indices
  hastet   = isfield(mesh, 'tet');  % tetraheders as a Mx4 matrix with vertex indices
  hashex   = isfield(mesh, 'hex');  % hexaheders  as a Mx8 matrix with vertex indices
end

% convert string into boolean values
faceindex   = istrue(faceindex);   % yes=view the face number
vertexindex = istrue(vertexindex); % yes=view the vertex number

if isempty(vertexcolor)
  if haspos && hascolor && (hastri || hastet || hashex || hasline || haspoly)
    vertexcolor = mesh.color;
  elseif haspos && (hastri || hastet || hashex || hasline || haspoly)
    vertexcolor ='none';
  else
    vertexcolor ='k';
  end
end

% there are various ways of specifying that this should not be plotted
if isequal(vertexcolor, 'false') || isequal(vertexcolor, 'no') || isequal(vertexcolor, 'off') || isequal(vertexcolor, false)
  vertexcolor = 'none';
end
if isequal(facecolor, 'false') || isequal(facecolor, 'no') || isequal(facecolor, 'off') || isequal(facecolor, false)
  facecolor = 'none';
end
if isequal(edgecolor, 'false') || isequal(edgecolor, 'no') || isequal(edgecolor, 'off') || isequal(edgecolor, false)
  edgecolor = 'none';
end

% color management
if ischar(vertexcolor) && exist([vertexcolor '.m'], 'file')
  vertexcolor = eval(vertexcolor);
elseif ischar(vertexcolor) && isequal(vertexcolor, 'curv') % default of ft_sourceplot method surface
  if isfield(mesh, 'curv')
    cortex_light = eval('cortex_light');
    cortex_dark  = eval('cortex_dark');
    % the curvature determines the color of gyri and sulci
    vertexcolor = mesh.curv(:) * cortex_dark + (1-mesh.curv(:)) * cortex_light;
  else
    cortex_light = eval('cortex_light');
    vertexcolor = repmat(cortex_light, size(mesh.pos,1), 1);
    warning('no curv field present in the mesh structure, using cortex_light as vertexcolor')
  end
end
if ischar(facecolor) && exist([facecolor '.m'], 'file')
  facecolor = eval(facecolor);
end
if ischar(edgecolor) && exist([edgecolor '.m'], 'file')
  edgecolor = eval(edgecolor);
end

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

if isfield(mesh, 'pos')
  % this is assumed to reflect 3-D vertices
  pos = mesh.pos;
elseif isfield(mesh, 'prj')
  % this happens sometimes if the 3-D vertices are projected to a 2-D plane
  pos = mesh.prj;
else
  error('no vertices found');
end

if isempty(pos)
  hs=[];
  return
end

if hastri+hastet+hashex+hasline+haspoly>1
  error('cannot deal with simultaneous triangles, tetraheders and/or hexaheders')
end

if hastri
  tri = mesh.tri;
elseif haspoly
  % these are treated just like triangles
  tri = mesh.poly;
elseif hastet
  % represent the tetraeders as the four triangles
  tri = [
    mesh.tet(:,[1 2 3]);
    mesh.tet(:,[2 3 4]);
    mesh.tet(:,[3 4 1]);
    mesh.tet(:,[4 1 2])];
  % or according to SimBio:  (1 2 3), (2 4 3), (4 1 3), (1 4 2)
  % there are shared triangles between neighbouring tetraeders, remove these
  tri = unique(tri, 'rows');
elseif hashex
  % represent the hexaheders as a collection of 6 patches
  tri = [
    mesh.hex(:,[1 2 3 4]);
    mesh.hex(:,[5 6 7 8]);
    mesh.hex(:,[1 2 6 5]);
    mesh.hex(:,[2 3 7 6]);
    mesh.hex(:,[3 4 8 7]);
    mesh.hex(:,[4 1 5 8]);
    ];
  % there are shared faces between neighbouring hexaheders, remove these
  tri = unique(tri, 'rows');
else
  tri = [];
end

if hasline
  line = mesh.line;
else
  line = [];
end

if haspos
  if ~isempty(tri)
    hs = patch('Vertices', pos, 'Faces', tri);
  elseif ~isempty(line)
    hs = patch('Vertices', pos, 'Faces', line);
  else
    hs = patch('Vertices', pos, 'Faces', []);
  end
  %set(hs, 'FaceColor', facecolor);
  set(hs, 'EdgeColor', edgecolor);
  set(hs, 'tag', tag);
end

% the vertexcolor can be specified either as a RGB color for each vertex, or as a single value at each vertex
% the facecolor can be specified either as a RGB color for each triangle, or as a single value at each triangle
% if there are triangles, the vertexcolor is used for linear interpolation over the patches
vertexpotential = ~isempty(tri) && ~ischar(vertexcolor) && (size(pos,1)==numel(vertexcolor) || size(pos,1)==size(vertexcolor,1) && (size(vertexcolor,2)==1 || size(vertexcolor,2)==3));
facepotential   = ~isempty(tri) && ~ischar(facecolor  ) && (size(tri,1)==numel(facecolor  ) || size(tri,1)==size(facecolor  ,1) && (size(facecolor  ,2)==1 || size(facecolor,  2)==3));

% if both vertexcolor and facecolor are numeric arrays, let the vertexcolor prevail
if vertexpotential
  % vertexcolor is an array with number of elements equal to the number of vertices
  set(hs, 'FaceVertexCData', vertexcolor, 'FaceColor', 'interp');
elseif facepotential
  set(hs, 'FaceVertexCData', facecolor, 'FaceColor', 'flat');
else
  % the color is indicated as a single character or as a single RGB triplet
  set(hs, 'FaceColor', facecolor);
end

% if facealpha is an array with number of elements equal to the number of vertices
if size(pos,1)==numel(facealpha)
  set(hs, 'FaceVertexAlphaData', facealpha);
  set(hs, 'FaceAlpha', 'interp');
elseif ~isempty(pos) && numel(facealpha)==1 && facealpha~=1
  % the default is 1, so that does not have to be set
  set(hs, 'FaceAlpha', facealpha);
end

if edgealpha~=1
  % the default is 1, so that does not have to be set
  set(hs, 'EdgeAlpha', edgealpha);
end

if faceindex
  % plot the triangle indices (numbers) at each face
  for face_indx=1:size(tri,1)
    str = sprintf('%d', face_indx);
    tri_x = (pos(tri(face_indx,1), 1) +  pos(tri(face_indx,2), 1) +  pos(tri(face_indx,3), 1))/3;
    tri_y = (pos(tri(face_indx,1), 2) +  pos(tri(face_indx,2), 2) +  pos(tri(face_indx,3), 2))/3;
    tri_z = (pos(tri(face_indx,1), 3) +  pos(tri(face_indx,2), 3) +  pos(tri(face_indx,3), 3))/3;
    h   = text(tri_x, tri_y, tri_z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hs  = [hs; h];
  end
end

if ~isequal(vertexcolor, 'none') && ~vertexpotential
  % plot the vertices as points
  
  if isempty(vertexcolor)
    % use black for all points
    if isscalar(vertexsize)
      if size(pos,2)==2
        hs = plot(pos(:,1), pos(:,2), ['k' vertexmarker]);
      else
        hs = plot3(pos(:,1), pos(:,2), pos(:,3), ['k' vertexmarker]);
      end
      set(hs, 'MarkerSize', vertexsize);
    else
      if size(pos,2)==2
        for i=1:size(pos,1)
          hs = plot(pos(i,1), pos(i,2), ['k' vertexmarker]);
          set(hs, 'MarkerSize', vertexsize(i));
        end
      else
        for i=1:size(pos,1)
          hs = plot3(pos(i,1), pos(i,2), pos(i,3), ['k' vertexmarker]);
          set(hs, 'MarkerSize', vertexsize(i));
        end
      end
    end
    
  elseif ischar(vertexcolor) && numel(vertexcolor)==1
    % one color for all points
    if isscalar(vertexsize)
      if size(pos,2)==2
        hs = plot(pos(:,1), pos(:,2), [vertexcolor vertexmarker]);
      else
        hs = plot3(pos(:,1), pos(:,2), pos(:,3), [vertexcolor vertexmarker]);
      end
      set(hs, 'MarkerSize', vertexsize);
    else
      if size(pos,2)==2
        for i=1:size(pos,1)
          hs = plot(pos(i,1), pos(i,2), [vertexcolor vertexmarker]);
          set(hs, 'MarkerSize', vertexsize(i));
        end
      else
        for i=1:size(pos,1)
          hs = plot3(pos(i,1), pos(i,2), pos(i,3), [vertexcolor vertexmarker]);
          set(hs, 'MarkerSize', vertexsize(i));
        end
      end
    end
    
  elseif ischar(vertexcolor) && numel(vertexcolor)==size(pos,1)
    % one color for each point
    if size(pos,2)==2
      for i=1:size(pos,1)
        hs = plot(pos(i,1), pos(i,2), [vertexcolor(i) vertexmarker]);
        if isscalar(vertexsize)
          set(hs, 'MarkerSize', vertexsize);
        else
          set(hs, 'MarkerSize', vertexsize(i));
        end
      end
    else
      for i=1:size(pos,1)
        hs = plot3(pos(i,1), pos(i,2), pos(i,3), [vertexcolor(i) vertexmarker]);
        if isscalar(vertexsize)
          set(hs, 'MarkerSize', vertexsize);
        else
          set(hs, 'MarkerSize', vertexsize(i));
        end
      end
    end
    
  elseif ~ischar(vertexcolor) && size(vertexcolor,1)==1
    % one RGB color for all points
    if size(pos,2)==2
      hs = plot(pos(:,1), pos(:,2), vertexmarker);
      set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor);
    else
      hs = plot3(pos(:,1), pos(:,2), pos(:,3), vertexmarker);
      set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor);
    end
    
  elseif ~ischar(vertexcolor) && size(vertexcolor,1)==size(pos,1) && size(vertexcolor,2)==3
    % one RGB color for each point
    if size(pos,2)==2
      for i=1:size(pos,1)
        hs = plot(pos(i,1), pos(i,2), vertexmarker);
        if isscalar(vertexsize)
          set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor(i,:));
        else
          set(hs, 'MarkerSize', vertexsize(i), 'MarkerEdgeColor', vertexcolor(i,:));
        end
      end
    else
      for i=1:size(pos,1)
        hs = plot3(pos(i,1), pos(i,2), pos(i,3), vertexmarker);
        if isscalar(vertexsize)
          set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor(i,:));
        else
          set(hs, 'MarkerSize', vertexsize(i), 'MarkerEdgeColor', vertexcolor(i,:));
        end
      end
    end
    
  else
    error('Unknown color specification for the vertices');
  end
  
end % plotting the vertices as points

if vertexindex
  % plot the vertex indices (numbers) at each node
  for node_indx=1:size(pos,1)
    str = sprintf('%d', node_indx);
    if size(pos, 2)==2
      h = text(pos(node_indx, 1), pos(node_indx, 2), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    else
      h = text(pos(node_indx, 1), pos(node_indx, 2), pos(node_indx, 3), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    hs = [hs; h];
  end
end

axis off
axis vis3d
axis equal

if ~nargout
  clear hs
end
if ~holdflag
  hold off
end

warning(ws); % revert to original state
