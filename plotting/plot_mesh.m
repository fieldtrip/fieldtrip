function plot_mesh(bnd, varargin)

% PLOT_MESH visualizes the information of a mesh contained in the first
% argument bnd. The boundary argument (bnd) contains typically 2 fields
% called .pnt and .tri referring to vertices and triangulation of a mesh.
%
% Use as
%   plot_mesh(bnd, ...)
%
% PLOT_MESH also allows to plot only vertices by
%   plot_mesh(pnt)
% where pnt is a list of 3d points cartesian coordinates.
%
% Graphic facilities are available for vertices, edges and faces. A list of
% the arguments is given below with the correspondent admitted choices.
%
%     'facecolor'     [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'vertexcolor'   [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'edgecolor'     [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'faceindex'     true or false
%     'vertexindex'   true or false
%     'facealpha'     transparency, between 0 and 1
%
% If you don't want the faces or vertices to be plotted, you should
% specify facecolor or respectively edgecolor as 'none'.
%
% Example
%   [pnt, tri] = icosahedron162;
%   bnd.pnt = pnt;
%   bnd.tri = tri;
%   plot_mesh(bnd, 'facecolor', 'skin', 'edgecolor', 'none')
%   camlight
%
% See also TRIMESH

% Copyright (C) 2009, Cristiano Micheli
%
% $Log: plot_mesh.m,v $
% Revision 1.29  2009/09/24 07:46:50  crimic
% added realistic colors
%
% Revision 1.28  2009/09/23 09:27:35  roboos
% removed accidental keyboard statement
%
% Revision 1.27  2009/09/23 09:26:55  roboos
% empty is not an appropriate facecolor
%
% Revision 1.26  2009/09/23 09:18:17  roboos
% made handling of facecolor=none more robust, idem for vertex and edgecolor
%
% Revision 1.25  2009/09/23 06:53:20  roboos
% updated documentation, some cleanup of the code and defaults (no functional changes)
%
% Revision 1.24  2009/09/23 06:41:34  roboos
% added "see also"
%
% Revision 1.23  2009/07/29 06:40:16  roboos
% allow empty pnt
%
% Revision 1.22  2009/06/29 16:03:16  roboos
% allow input Nx3 as set of points
%
% Revision 1.21  2009/06/25 16:02:03  crimic
% fixed little error
%
% Revision 1.20  2009/06/21 19:45:52  crimic
% minor changes
%
% Revision 1.19  2009/06/21 19:25:01  crimic
% added tag argument and set edge default color to black
%
% Revision 1.18  2009/06/16 12:19:04  crimic
% erased output graphic handles
%
% Revision 1.17  2009/06/15 15:48:21  roboos
% fixed handling in case input does not have a triangulation
%
% Revision 1.16  2009/06/08 11:52:55  crimic
% updates single points' vertexcolor property
%
% Revision 1.15  2009/06/04 12:55:34  crimic
% changes input parameters check
%
% Revision 1.14  2009/06/03 09:53:50  roboos
% added option for facealpha (transparency)
% updated help for consistency
%
% Revision 1.13  2009/05/13 07:54:36  crimic
% updated help
%
% Revision 1.12  2009/05/13 07:49:55  crimic
% inserted option to manage plot of points in 3d
%
% Revision 1.11  2009/05/12 18:11:21  roboos
% cleaned up whitespace and indentation
%
% Revision 1.10  2009/04/17 13:43:33  crimic
% updated help
%
% Revision 1.9  2009/04/17 13:38:23  crimic
% added default options for varargin, added edgecolor argument
%
% Revision 1.8  2009/04/09 09:36:50  crimic
% *** empty log message ***
%
% Revision 1.7  2009/04/09 09:35:53  crimic
% integrated help
%
% Revision 1.6  2009/04/09 09:20:33  crimic
% clean-up of non used options (val, contour and surface), inserted istrue check
%
% Revision 1.5  2009/04/08 17:09:40  crimic
% added help and modified input argument structure
%
% Revision 1.4  2009/04/08 16:08:50  crimic
% indented with 2 spaces tab, added log signature and copyright
%

% FIXME: introduce option for color coding (see sourceplot)
keyvalcheck(varargin, 'forbidden', {'faces', 'edges', 'vertices'});

if ~isstruct(bnd) && isnumeric(bnd) && size(bnd,2)==3
  % the input seems like a list of points, convert into something that resembles a mesh
  warning('off', 'MATLAB:warn_r14_stucture_assignment');
  bnd.pnt = bnd;
end

% get the optional input arguments
facecolor   = keyval('facecolor',   varargin); if isempty(facecolor),   facecolor='white';end
vertexcolor = keyval('vertexcolor', varargin); if isempty(vertexcolor), vertexcolor='none';end
edgecolor   = keyval('edgecolor',   varargin); if isempty(edgecolor),   edgecolor='k';end
faceindex   = keyval('faceindex',   varargin); if isempty(faceindex),   faceindex=false;end
vertexindex = keyval('vertexindex', varargin); if isempty(vertexindex), vertexindex=false;end
vertexsize  = keyval('vertexsize',  varargin); if isempty(vertexsize),  vertexsize=10;end
facealpha   = keyval('facealpha',   varargin); if isempty(facealpha),   facealpha=1;end
tag         = keyval('tag',         varargin); if isempty(tag),         tag='';end

% convert string into boolean values
faceindex   = istrue(faceindex);
vertexindex = istrue(vertexindex);

skin   = [255 213 119]/255;
brain  = [202 100 100]/255;
cortex = [255 213 119]/255;
    
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
hold on

if ~isfield(bnd, 'tri')
  bnd.tri = [];
end

pnt = bnd.pnt;
tri = bnd.tri;

if ~isempty(pnt)
  hs = patch('Vertices', pnt, 'Faces', tri);
  set(hs, 'FaceColor', facecolor);
  set(hs, 'FaceAlpha', facealpha);
  set(hs, 'EdgeColor', edgecolor);
  set(hs, 'tag', tag);
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

if ~isequal(vertexcolor, 'none')
  if size(pnt, 2)==2
    hs = plot(pnt(:,1), pnt(:,2), 'k.');
  else
    hs = plot3(pnt(:,1), pnt(:,2), pnt(:,3), 'k.');
  end
  if ~isempty(vertexcolor)
    try
      set(hs, 'Marker','.','MarkerEdgeColor', vertexcolor,'MarkerSize', vertexsize);
    catch
      error('Unknown color')
    end
  end
  if vertexindex
    % plot the vertex indices (numbers) at each node
    for node_indx=1:size(pnt,1)
      str = sprintf('%d', node_indx);
      if size(pnt, 2)==2
        h   = text(pnt(node_indx, 1), pnt(node_indx, 2), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
      else
        h   = text(pnt(node_indx, 1), pnt(node_indx, 2), pnt(node_indx, 3), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
      end
      hs  = [hs; h];
    end
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

