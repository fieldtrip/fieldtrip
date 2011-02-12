function [hs, hc, contour] = triplot(pnt, tri, val, mode, levels)

% TRIPLOT make 2D or 3D plot of triangulated surface and interpolated values
% The surface can be displayed with linear interpolated values or with
% linear interpolated contours of a potential distribution.
%
% Use as
%   triplot(pnt, tri, value)
%   triplot(pnt, tri, value, mode)
%   triplot(pnt, tri, value, mode, levels)
% This will make a plot of value on the surface described by triangles
% tri with vertices pnt. The matrix tri can be [], in which case a
% it will be computed using a delaunay triangulation.
%
% The visualization mode can be
%   'surface'     make interpolated plot of value on surface (default)
%   'faces'       plot white triangles only        (value can be [])
%   'faces_red'   plot red triangles only          (value can be [])
%   'faces_blue'  plot blue triangles only         (value can be [])
%   'faces_skin'  plot skin-colored triangles only (value can be [])
%   'face_index'  plot index of each triangle      (value can be [])
%   'nodes'       plot black vertices only         (value can be [])
%   'node_index'  plot index of each vertex        (value can be [])
%   'node_label'  plot label of each vertex        (value should be cell array)
%   'edges'       plot black edges only            (value can be [])
%   'contour'     make interpolated contour plot of value on surface
%   'contour_bw'  make interpolated black-white contour plot
%   'contour_rb'  make interpolated contour plot with red-blue
%
% With the optional levels, you can specify the levels at which contours will
% be plotted
%
% See also PATCH, COLORMAP, VIEW (general Matlab commands)

% updated on Mon Jul 23 12:41:44 MET DST 2001
% updated on Fri Jan 31 11:47:28 CET     2003

% Copyright (C) 2001=2006, Robert Oostenveld
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
% $Id: triplot.m 952 2010-04-21 18:29:51Z roboos $

% start with empty return values
hs      = [];
hc      = [];
contour = [];

% everything is added to the current figure
holdflag = ishold;
hold on

% check the input variables
if ~isempty(val)
  val = val(:);
end

if nargin<4
  mode = 'surface';
end

if isempty(tri) & ~strcmp(mode, 'nodes')
  % no triangulation was specified but a triangulation is needed
  if size(pnt,2)==2
    % make a 2d triangulation of the points using delaunay
    tri = delaunay(pnt(:,1), pnt(:,2));
  else
    % make a 2d triangulation of the projected points using delaunay
    prj = elproj(pnt);
    tri = delaunay(prj(:,1), prj(:,2));
  end
end

if size(tri,2)==2
  % lines are specified instead of triangles, convert to triangles
  tri(:,3) = tri(:,2);
end

if nargin<3
  warning('only displaying triangle edges')
  mode='edges';
  val = [];
elseif nargin<4
  % warning('default displaying surface')
  mode='surface';
elseif nargin<5
  % determine contour levels
  if ~isempty(val) & ~iscell(val)
    absmax = max(abs([min(val) max(val)]));
    levels = linspace(-absmax,absmax,21);
  else
    levels = [];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute contours for 2D or 3D triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'contour') || strcmp(mode, 'contour_bw') || strcmp(mode, 'contour_rb')
  triangle_val = val(tri);
  triangle_min = min(triangle_val, [], 2);
  triangle_max = max(triangle_val, [], 2);

  for cnt_indx=1:length(levels)
    cnt = levels(cnt_indx);
    use = cnt>=triangle_min & cnt<=triangle_max;
    counter = 0;
    intersect1 = [];
    intersect2 = [];

    for tri_indx=find(use)'
      pos  = pnt(tri(tri_indx,:), :);
      v(1) = triangle_val(tri_indx,1);
      v(2) = triangle_val(tri_indx,2);
      v(3) = triangle_val(tri_indx,3);
      la(1) = (cnt-v(1)) / (v(2)-v(1)); % abcissa between vertex 1 and 2
      la(2) = (cnt-v(2)) / (v(3)-v(2)); % abcissa between vertex 2 and 3
      la(3) = (cnt-v(3)) / (v(1)-v(3)); % abcissa between vertex 1 and 2
      abc(1,:) = pos(1,:) + la(1) * (pos(2,:) - pos(1,:));
      abc(2,:) = pos(2,:) + la(2) * (pos(3,:) - pos(2,:));
      abc(3,:) = pos(3,:) + la(3) * (pos(1,:) - pos(3,:));
      counter = counter + 1;
      sel     = find(la>=0 & la<=1);
      intersect1(counter, :) = abc(sel(1),:);
      intersect2(counter, :) = abc(sel(2),:);
    end

    % remember the details for external reference
    contour(cnt_indx).level = cnt;
    contour(cnt_indx).n     = counter;
    contour(cnt_indx).intersect1 = intersect1;
    contour(cnt_indx).intersect2 = intersect2;
  end

  % collect all different contourlevels for plotting
  intersect1 = [];
  intersect2 = [];
  cntlevel   = [];
  for cnt_indx=1:length(levels)
    intersect1 = [intersect1; contour(cnt_indx).intersect1];
    intersect2 = [intersect2; contour(cnt_indx).intersect2];
    cntlevel   = [cntlevel; ones(contour(cnt_indx).n,1) * levels(cnt_indx)];
  end

  X = [intersect1(:,1) intersect2(:,1)]';
  Y = [intersect1(:,2) intersect2(:,2)]';
  C = [cntlevel(:)     cntlevel(:)]';

  if size(pnt,2)>2
    Z = [intersect1(:,3) intersect2(:,3)]';
  else
    Z = zeros(2, length(cntlevel));
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the desired detail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch lower(mode)

  case 'faces'
    % plot the faces of the 2D or 3D triangulation
    hs = patch('Vertices', pnt, 'Faces', tri);
    set(hs, 'FaceColor', 'white');
    set(hs, 'EdgeColor', 'none');

  case 'faces_skin'
    % plot the faces of the 2D or 3D triangulation
    skin   = [255 213 119]/255;
    brain  = [202 100 100]/255;
    cortex = [255 213 119]/255;
    hs = patch('Vertices', pnt, 'Faces', tri);
    set(hs, 'FaceColor', skin);
    set(hs, 'EdgeColor', 'none');
    lighting gouraud
    material shiny
    camlight

  case 'faces_red'
    % plot the faces of the 2D or 3D triangulation
    hs = patch('Vertices', pnt, 'Faces', tri);
    set(hs, 'FaceColor', [1 0 0]);
    set(hs, 'EdgeColor', 'none');
    lighting gouraud
    material shiny
    camlight
    
  case 'faces_blue'
    % plot the faces of the 2D or 3D triangulation
    hs = patch('Vertices', pnt, 'Faces', tri);
    set(hs, 'FaceColor', [0 0 1]);
    set(hs, 'EdgeColor', 'none');
    lighting gouraud
    material shiny
    camlight

  case 'face_index'
    % plot the triangle indices (numbers) at each face
    for face_indx=1:size(tri,1)
      str = sprintf('%d', face_indx);
      tri_x = (pnt(tri(face_indx,1), 1) +  pnt(tri(face_indx,2), 1) +  pnt(tri(face_indx,3), 1))/3;
      tri_y = (pnt(tri(face_indx,1), 2) +  pnt(tri(face_indx,2), 2) +  pnt(tri(face_indx,3), 2))/3;
      tri_z = (pnt(tri(face_indx,1), 3) +  pnt(tri(face_indx,2), 3) +  pnt(tri(face_indx,3), 3))/3;
      h   = text(tri_x, tri_y, tri_z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
      hs  = [hs; h];
    end

  case 'edges'
    % plot the edges of the 2D or 3D triangulation
    hs = patch('Vertices', pnt, 'Faces', tri);
    set(hs, 'FaceColor', 'none');
    set(hs, 'EdgeColor', 'black');

  case 'nodes'
    % plot the nodes (vertices) only as points
    if size(pnt, 2)==2
      hs = plot(pnt(:,1), pnt(:,2), 'k.');
    else
      hs = plot3(pnt(:,1), pnt(:,2), pnt(:,3), 'k.');
    end
    
  case 'nodes_blue'
    % plot the nodes (vertices) only as points
    if size(pnt, 2)==2
      hs = plot(pnt(:,1), pnt(:,2), 'b.', 'MarkerSize', 20);
    else
      hs = plot3(pnt(:,1), pnt(:,2), pnt(:,3), 'b.', 'MarkerSize', 20);
    end  
    
  case 'nodes_red'
    % plot the nodes (vertices) only as points
    if size(pnt, 2)==2
      hs = plot(pnt(:,1), pnt(:,2), 'r.', 'MarkerSize', 20);
    else
      hs = plot3(pnt(:,1), pnt(:,2), pnt(:,3), 'r.', 'MarkerSize', 20);
    end   

  case 'node_index'
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

  case 'node_label'
    % plot the vertex indices (numbers) at each node
    for node_indx=1:size(pnt,1)
      str = val{node_indx};
      if ~isempty(str)
        if size(pnt, 2)==2
          h   = text(pnt(node_indx, 1), pnt(node_indx, 2), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        else
          h   = text(pnt(node_indx, 1), pnt(node_indx, 2), pnt(node_indx, 3), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
      else
        h = -1;
      end
      hs  = [hs; h];
    end

  case 'surface'
    % plot a 2D or 3D triangulated surface with linear interpolation
    if length(val)==size(pnt,1)
      hs = patch('Vertices', pnt, 'Faces', tri, 'FaceVertexCData', val, 'FaceColor', 'interp');
    else
      hs = patch('Vertices', pnt, 'Faces', tri, 'CData', val, 'FaceColor', 'flat');
    end
    set(hs, 'EdgeColor', 'none');

  case 'contour_bw'
    % make black-white contours
    hc = [];
    for i=1:length(cntlevel)
      if cntlevel(i)>0
        linestyle = '-';
        linewidth = 1;
      elseif cntlevel(i)<0
        linestyle = '--';
        linewidth = 1;
      else
        linestyle = '-';
        linewidth = 2;
      end
      h1 = patch('XData', X(:,i), 'Ydata', Y(:,i), ...
        'ZData', Z(:,i), 'CData', C(:,i), ...
        'facecolor','none','edgecolor','black', ...
        'linestyle', linestyle, 'linewidth', linewidth, ...
        'userdata',cntlevel(i));
      hc = [hc; h1];
    end

  case 'contour_rb'
    % make red-blue contours
    hc = [];
    for i=1:length(cntlevel)
      if cntlevel(i)>0
        edgecolor = 'red';
      elseif cntlevel(i)<0
        edgecolor = 'blue';
      else
        edgecolor = 'black';
      end
      h1 = patch('XData', X(:,i), 'Ydata', Y(:,i), ...
        'ZData', Z(:,i), 'CData', C(:,i), ...
        'facecolor','none','edgecolor',edgecolor, ...
        'linestyle', '-', 'linewidth', 3, ...
        'userdata',cntlevel(i));
      hc = [hc; h1];
    end

  case 'contour'
    % make full-color contours
    hc = [];
    for i=1:length(cntlevel)
      h1 = patch('XData', X(:,i), 'Ydata', Y(:,i), ...
        'ZData', Z(:,i), 'CData', C(:,i), ...
        'facecolor','none','edgecolor','flat',...
        'userdata',cntlevel(i));
      hc = [hc; h1];
    end

end % switch

axis off
axis vis3d
axis equal

if nargout==0
  clear contour hc hs
end

if ~holdflag
  hold off
end

