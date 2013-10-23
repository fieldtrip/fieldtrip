function [hs] = ft_plot_mesh(bnd, varargin)

% FT_PLOT_MESH visualizes the information of a mesh contained in the first
% argument bnd. The boundary argument (bnd) typically contains two fields
% called .pnt and .tri referring to the vertices and the triangulation of
% the mesh.
%
% Use as
%   ft_plot_mesh(bnd, ...)
% or if you only want to plot the 3-D vertices
%   ft_plot_mesh(pnt, ...)
%
% Optional arguments should come in key-value pairs and can include
%     'facecolor'   = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'vertexcolor' = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r', or an Nx1 array where N is the number of vertices
%     'edgecolor'   = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'faceindex'   = true or false
%     'vertexindex' = true or false
%     'facealpha'   = transparency, between 0 and 1 (default = 1)
%     'edgealpha'   = transparency, between 0 and 1 (default = 1)
%     'surfaceonly' = true or false, plot only the outer surface of a hexahedral or tetrahedral mesh (default = false)
%
% If you don't want the faces or vertices to be plotted, you should specify the color as 'none'.
%
% Example
%   [pnt, tri] = icosahedron162;
%   bnd.pnt = pnt;
%   bnd.tri = tri;
%   ft_plot_mesh(bnd, 'facecolor', 'skin', 'edgecolor', 'none')
%   camlight
%
% See also TRIMESH, PATCH

% Copyright (C) 2009-2013, Robert Oostenveld
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

if ~isstruct(bnd) && isnumeric(bnd) && size(bnd,2)==3
  % the input seems like a list of points, convert into something that resembles a mesh
  warning('off', 'MATLAB:warn_r14_stucture_assignment');
  bnd.pnt = bnd;
elseif isfield(bnd, 'pos')
  % the input seems to be a set of points from ft_prepare_sourcemodel or ft_dipolefitting
  bnd.pnt = bnd.pos;
end

% the input is a structure, but might also be a struct-array
if numel(bnd)>1
  % plot each of the boundaries
  for i=1:numel(bnd)
    ft_plot_mesh(bnd(i), varargin{:})
  end
  return
end

% get the optional input arguments
vertexcolor = ft_getopt(varargin, 'vertexcolor');
facecolor   = ft_getopt(varargin, 'facecolor',   'white');
edgecolor   = ft_getopt(varargin, 'edgecolor',   'k');
faceindex   = ft_getopt(varargin, 'faceindex',   false);
vertexindex = ft_getopt(varargin, 'vertexindex', false);
vertexsize  = ft_getopt(varargin, 'vertexsize',  10);
facealpha   = ft_getopt(varargin, 'facealpha',   1);
edgealpha   = ft_getopt(varargin, 'edgealpha',   1);
tag         = ft_getopt(varargin, 'tag',         '');
surfaceonly    = ft_getopt(varargin, 'surfaceonly',  false);

if surfaceonly
  bnd = mesh2edge(bnd);
end

haspnt  = isfield(bnd, 'pnt');  % vertices
hastri  = isfield(bnd, 'tri');  % triangles   as a Mx3 matrix with vertex indices
hastet  = isfield(bnd, 'tet');  % tetraheders as a Mx4 matrix with vertex indices
hashex  = isfield(bnd, 'hex');  % hexaheders  as a Mx8 matrix with vertex indices
hasline = isfield(bnd, 'line'); % line segments in 3-D
haspoly = isfield(bnd, 'poly'); % polynomial surfaces in 3-D

if (hastet || hashex) && ~surfaceonly
  warning('you probably want to use the "surfaceonly" option for plotting only the outer surface')
end

if isempty(vertexcolor)
  if haspnt && (hastri || hastet || hashex || hasline || haspoly)
    vertexcolor ='none';
  else
    vertexcolor ='k';
  end
end

% convert string into boolean values
faceindex   = istrue(faceindex);   % yes=view the face number
vertexindex = istrue(vertexindex); % yes=view the vertex number

% there a various ways of disabling the plotting
if isequal(vertexcolor, 'false') || isequal(vertexcolor, 'no') || isequal(vertexcolor, 'off') || isequal(vertexcolor, false)
  vertexcolor = 'none';
end
if isequal(facecolor, 'false') || isequal(facecolor, 'no') || isequal(facecolor, 'off') || isequal(facecolor, false)
  facecolor = 'none';
end
if isequal(edgecolor, 'false') || isequal(edgecolor, 'no') || isequal(edgecolor, 'off') || isequal(edgecolor, false)
  edgecolor = 'none';
end

% new colors management
if strcmpi(vertexcolor,'skin') || strcmpi(vertexcolor,'brain') || strcmpi(vertexcolor,'cortex')
  vertexcolor = eval(vertexcolor);
end
if strcmpi(facecolor,'skin') || strcmpi(facecolor,'brain') || strcmpi(facecolor,'cortex')
  facecolor = eval(facecolor);
end

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

if isfield(bnd, 'pnt')
  % this is assumed to reflect 3-D vertices
  pnt = bnd.pnt;
elseif isfield(bnd, 'prj')
  % this happens sometimes if the 3-D vertices are projected to a 2-D plane
  pnt = bnd.prj;
else
  error('no vertices found');
end

if isempty(pnt)
    hs=[];
    return
end

if hastri+hastet+hashex+hasline+haspoly>1
  error('cannot deal with simultaneous triangles, tetraheders and/or hexaheders')
end

if hastri
  tri = bnd.tri;
elseif haspoly
  % these are treated just like triangles
  tri = bnd.poly;
elseif hastet
  % represent the tetraeders as the four triangles
  tri = [
    bnd.tet(:,[1 2 3]);
    bnd.tet(:,[2 3 4]);
    bnd.tet(:,[3 4 1]);
    bnd.tet(:,[4 1 2])];
  % or according to SimBio:  (1 2 3), (2 4 3), (4 1 3), (1 4 2)
  % there are shared triangles between neighbouring tetraeders, remove these
  tri = unique(tri, 'rows');
elseif hashex
  % represent the hexaheders as a collection of 6 patches
  tri = [
    bnd.hex(:,[1 2 3 4]);
    bnd.hex(:,[5 6 7 8]);
    bnd.hex(:,[1 2 6 5]);
    bnd.hex(:,[2 3 7 6]);
    bnd.hex(:,[3 4 8 7]);
    bnd.hex(:,[4 1 5 8]);
    ];
  % there are shared faces between neighbouring hexaheders, remove these
  tri = unique(tri, 'rows');
else
  tri = [];
end

if hasline
  line = bnd.line;
else
  line = [];
end

if haspnt 
  if ~isempty(tri)
    hs = patch('Vertices', pnt, 'Faces', tri);
  elseif ~isempty(line)
    hs = patch('Vertices', pnt, 'Faces', line);
  else
    hs = patch('Vertices', pnt, 'Faces', []);
  end
  %set(hs, 'FaceColor', facecolor);
  set(hs, 'EdgeColor', edgecolor);
  set(hs, 'tag', tag);
end

% the vertexcolor can be specified either as a color for each point that will be drawn, or as a value at each vertex
% if there are triangles, the vertexcolor is used for linear interpolation over the patches
vertexpotential = ~isempty(tri) && ~ischar(vertexcolor) && (size(pnt,1)==numel(vertexcolor) || size(pnt,1)==size(vertexcolor,1));
facepotential   = ~isempty(tri) && ~ischar(facecolor)   && (size(tri,1)==numel(facecolor)   || size(tri,1)==size(facecolor,1));
if facepotential
  set(hs, 'FaceVertexCData', facecolor, 'FaceColor', 'flat');
else
  set(hs, 'FaceColor', facecolor);
end

if vertexpotential
  % vertexcolor is an array with number of elements equal to the number of vertices
  % if both vertexcolor and facecolor are arrays, let the vertexcolor
  % prevail
  set(hs, 'FaceVertexCData', vertexcolor, 'FaceColor', 'interp');
end

% if facealpha is an array with number of elements equal to the number of vertices
if size(pnt,1)==numel(facealpha)
  set(hs, 'FaceVertexAlphaData', facealpha);
  set(hs, 'FaceAlpha', 'interp');
elseif ~isempty(pnt) && numel(facealpha)==1 && facealpha~=1
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
    tri_x = (pnt(tri(face_indx,1), 1) +  pnt(tri(face_indx,2), 1) +  pnt(tri(face_indx,3), 1))/3;
    tri_y = (pnt(tri(face_indx,1), 2) +  pnt(tri(face_indx,2), 2) +  pnt(tri(face_indx,3), 2))/3;
    tri_z = (pnt(tri(face_indx,1), 3) +  pnt(tri(face_indx,2), 3) +  pnt(tri(face_indx,3), 3))/3;
    h   = text(tri_x, tri_y, tri_z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hs  = [hs; h];
  end
end

if ~isequal(vertexcolor, 'none') && ~vertexpotential
  % plot the vertices as points
  
  if isempty(vertexcolor)
    % use black for all points
    if size(pnt,2)==2
      hs = plot(pnt(:,1), pnt(:,2), 'k.');
    else
      hs = plot3(pnt(:,1), pnt(:,2), pnt(:,3), 'k.');
    end
    set(hs, 'MarkerSize', vertexsize);
    
  elseif ischar(vertexcolor) && numel(vertexcolor)==1
    % one color for all points
    if size(pnt,2)==2
      hs = plot(pnt(:,1), pnt(:,2), [vertexcolor '.']);
    else
      hs = plot3(pnt(:,1), pnt(:,2), pnt(:,3), [vertexcolor '.']);
    end
    set(hs, 'MarkerSize', vertexsize);
    
  elseif ischar(vertexcolor) && numel(vertexcolor)==size(pnt,1)
    % one color for each point
    if size(pnt,2)==2
      for i=1:size(pnt,1)
        hs = plot(pnt(i,1), pnt(i,2), [vertexcolor(i) '.']);
        set(hs, 'MarkerSize', vertexsize);
      end
    else
      for i=1:size(pnt,1)
        hs = plot3(pnt(i,1), pnt(i,2), pnt(i,3), [vertexcolor(i) '.']);
        set(hs, 'MarkerSize', vertexsize);
      end
    end
    
  elseif ~ischar(vertexcolor) && size(vertexcolor,1)==1
    % one RGB color for all points
    if size(pnt,2)==2
      hs = plot(pnt(:,1), pnt(:,2), '.');
      set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor);
    else
      hs = plot3(pnt(:,1), pnt(:,2), pnt(:,3), '.');
      set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor);
    end
    
  elseif ~ischar(vertexcolor) && size(vertexcolor,1)==size(pnt,1)
    % one RGB color for each point
    if size(pnt,2)==2
      for i=1:size(pnt,1)
        hs = plot(pnt(i,1), pnt(i,2), '.');
        set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor(i,:));
      end
    else
      for i=1:size(pnt,1)
        hs = plot3(pnt(i,1), pnt(i,2), pnt(i,3), '.');
        set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor(i,:));
      end
    end
    
  else
    error('Unknown color specification for the vertices');
  end
  
end % plotting the vertices as points

if vertexindex
  % plot the vertex indices (numbers) at each node
  for node_indx=1:size(pnt,1)
    str = sprintf('%d', node_indx);
    if size(pnt, 2)==2
      h = text(pnt(node_indx, 1), pnt(node_indx, 2), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    else
      h = text(pnt(node_indx, 1), pnt(node_indx, 2), pnt(node_indx, 3), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
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
