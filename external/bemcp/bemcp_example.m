
function bemcp_example
% Simple function to test/demonstrate how the Boundary element functions are
% used in combination with Fildtrip/Forwinv routines.
%
% 1. A model is created as 3 concentric meshed spheres (using FT's 
% icosahedron routines), 
% 2. then random electrodes are placed on the upper part of the outer 
% sphere. 
% 3. the model is then "prepared" with 'ft_prepare_bemmodel', this bits
% takes most time as it requires LOTS of calculation.
% 4. sensors and volumes are plugged together by 'forwinv_prepare_vol_sens'
% 5. Finally the leadfiled for 3 orthogonal sources placed at one location
% is calculated with 'forwinv_compute_leadfield.m'
% 6. Display the 3 leadfields
%
% NOTE: 
% this bit of code needs access to low level fieldtrip/forwinv routines 
% which have been copy/pasted here under.
% Be aware that this way of programming is generally NOT advisable!
% I used it only to ensure a quick & dirty check of the BEM module...

% Christophe Phillips
% $Id$

% create volume conductor starting from unit sphere
[pnt, tri] = mesh_sphere(162);

vol = [];
vol.cond = [1 1/80 1];
vol.source = 1; % index of source compartment
vol.skin_surface   = 3; % index of skin surface
% inner_skull_surface
vol.bnd(1).pnt = pnt*88;
vol.bnd(1).tri = tri;
% outer_skull_surface
vol.bnd(2).pnt = pnt*92;
vol.bnd(2).tri = tri;
% skin_surface
vol.bnd(3).pnt = pnt*100;
vol.bnd(3).tri = tri;

% create the BEM system matrix
cfg = [];
cfg.method = 'bemcp';
vol1 = ft_prepare_bemmodel(cfg, vol);


% create some random electrodes
pnt = randn(200,3);
pnt = pnt(pnt(:,3)>0, :);  % only those on the upper half
sens = [];
for i=1:size(pnt,1)
  sens.pnt(i,:) = pnt(i,:) / norm(pnt(i,:)); % scale towards the skin surface
  sens.label{i} = sprintf('%02d', i);
end

% prepare the sensor array and volume conduction, i.e. set up the linear
% interpolation from vertices to electrodes
[vol2, sens] = forwinv_prepare_vol_sens(vol1, sens);

lf = forwinv_compute_leadfield([0 0 50], sens, vol2);

figure; triplot(sens.pnt, [], lf(:,1)); colorbar
figure; triplot(sens.pnt, [], lf(:,2)); colorbar
figure; triplot(sens.pnt, [], lf(:,3)); colorbar

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions from FieldTrip
function [pnt, tri] = mesh_sphere(162)

% ICOSAHEDRON162 creates a 2-fold refined icosahedron

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: mesh_sphere(162).m,v $
% Revision 1.3  2003/11/28 09:40:12  roberto
% added a single help line
%
% Revision 1.2  2003/03/04 21:46:18  roberto
% added CVS log entry and synchronized all copyright labels
%

[pnt, tri] = icosahedron;
[pnt, tri] = refine(pnt, tri);
[pnt, tri] = refine(pnt, tri);

pnt = pnt ./ repmat(sqrt(sum(pnt.^2,2)), 1,3);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pnt, tri] = icosahedron

% ICOSAHEDRON creates an icosahedron
%
% [pnt, tri] = icosahedron
% creates an icosahedron with 12 vertices and 20 triangles
% 
% See also OCTAHEDRON, ICOSAHEDRON42, ICOSAHEDRON162, ICOSAHEDRON642, ICOSAHEDRON2562

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: icosahedron.m,v $
% Revision 1.4  2006/07/26 11:03:38  roboos
% added "see also octahedron"
%
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:18  roberto
% added CVS log entry and synchronized all copyright labels
%

tri = [
   1   2   3
   1   3   4
   1   4   5
   1   5   6
   1   6   2
   2   8   3
   3   9   4
   4  10   5
   5  11   6
   6   7   2
   7   8   2  
   8   9   3  
   9  10   4  
  10  11   5  
  11   7   6  
  12   8   7
  12   9   8
  12  10   9
  12  11  10
  12   7  11
];

pnt = zeros(12, 3);

rho=0.4*sqrt(5);
phi=2*pi*(0:4)/5;

pnt( 1, :) = [0 0  1];			% top point

pnt(2:6, 1) = rho*cos(phi)';
pnt(2:6, 2) = rho*sin(phi)';
pnt(2:6, 3) = rho/2;

pnt(7:11, 1) = rho*cos(phi - pi/5)';
pnt(7:11, 2) = rho*sin(phi - pi/5)';
pnt(7:11, 3) = -rho/2;

pnt(12, :) = [0 0 -1];			% bottom point

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pntr, tri] = refine(pnt, tri, method, varargin)

% REFINE a 3D surface that is described by a triangulation
%
% Use as
%   [pnt, tri] = refine(pnt, tri)
%   [pnt, tri] = refine(pnt, tri, 'updown', numtri)
%
% The default method is to refine the mesh globally by inserting a vertex at
% each edge according to the algorithm described in Banks, 1983.
%
% The alternative 'updown' method refines the mesh a couple of times
% using Banks' algorithm, followed by a downsampling using the REDUCEPATCH
% function.

% The Banks method is a memory efficient implementation which remembers
% the previously inserted vertices. The refinement algorithm executes in
% linear time with the number of triangles.

% Copyright (C) 2002-2005, Robert Oostenveld
%
% $Log: refine.m,v $
% Revision 1.4  2005/05/06 08:54:31  roboos
% added 'updown method, using matlab reducepatch function
%
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

if nargin<3
  method = 'banks';
end

switch lower(method)
case 'banks'
  npnt   = size(pnt,1);
  ntri   = size(tri,1);
  insert = spalloc(3*npnt,3*npnt,3*ntri);

  tri  = zeros(4*ntri,3);		% allocate memory for the new triangles
  pntr  = zeros(npnt+3*ntri,3);		% allocate memory for the maximum number of new vertices
  pntr(1:npnt,:) = pnt;			% insert the original vertices
  current = npnt;

  for i=1:ntri

    if ~insert(tri(i,1),tri(i,2))
      current = current + 1;
      pntr(current,:) = (pnt(tri(i,1),:) + pnt(tri(i,2),:))/2;
      insert(tri(i,1),tri(i,2)) = current;
      insert(tri(i,2),tri(i,1)) = current;
      v12 = current;
    else
      v12 = insert(tri(i,1),tri(i,2));
    end

    if ~insert(tri(i,2),tri(i,3))
      current = current + 1;
      pntr(current,:) = (pnt(tri(i,2),:) + pnt(tri(i,3),:))/2;
      insert(tri(i,2),tri(i,3)) = current;
      insert(tri(i,3),tri(i,2)) = current;
      v23 = current;
    else
      v23 = insert(tri(i,2),tri(i,3));
    end

    if ~insert(tri(i,3),tri(i,1))
      current = current + 1;
      pntr(current,:) = (pnt(tri(i,3),:) + pnt(tri(i,1),:))/2;
      insert(tri(i,3),tri(i,1)) = current;
      insert(tri(i,1),tri(i,3)) = current;
      v31 = current;
    else
      v31 = insert(tri(i,3),tri(i,1));
    end

    % add the 4 new triangles with the correct indices
    tri(4*(i-1)+1, :) = [tri(i,1) v12 v31];
    tri(4*(i-1)+2, :) = [tri(i,2) v23 v12];
    tri(4*(i-1)+3, :) = [tri(i,3) v31 v23];
    tri(4*(i-1)+4, :) = [v12 v23 v31];

  end

  % remove the space for the vertices that was not used
  pntr = pntr(1:current, :);

case 'updown'
  ntri = size(tri,1);
  while ntri<varargin{1}
    % increase the number of triangles by a factor of 4
    [pnt, tri] = refine(pnt, tri, 'banks');
    ntri = size(tri,1);
  end
  % reduce number of triangles using Matlab function
  [tri, pntr] = reducepatch(tri, pnt, varargin{1});

otherwise
  error(['unsupported method: ' method]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% $Log: triplot.m,v $
% Revision 1.7  2008/06/24 13:37:51  roboos
% added option faces_blue
% always return handle to objects that were plotted
%
% Revision 1.6  2007/01/03 17:00:35  roboos
% updated documentation, changed layout of code and comments
%
% Revision 1.5  2006/09/19 16:11:35  roboos
% added support for line segments, removed "axis equal" at end
%
% Revision 1.4  2006/05/02 19:15:28  roboos
% added 'faces_red' style
%
% Revision 1.3  2004/06/28 07:51:39  roberto
% improved documentation, added faces_skin
%
% Revision 1.2  2003/03/17 10:37:29  roberto
% improved general help comments and added copyrights
%

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
      la(1) = (cnt-v(1)) / (v(2)-v(1));	% abcissa between vertex 1 and 2
      la(2) = (cnt-v(2)) / (v(3)-v(2));	% abcissa between vertex 2 and 3
      la(3) = (cnt-v(3)) / (v(1)-v(3));	% abcissa between vertex 1 and 2
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
    skin_surface   = [255 213 119]/255;
    inner_skull_surface  = [202 100 100]/255;
    cortex = [255 213 119]/255;
    hs = patch('Vertices', pnt, 'Faces', tri);
    set(hs, 'FaceColor', skin_surface);
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

end	% switch

axis off
axis vis3d
axis equal

if nargout==0
  clear contour hc hs
end

if ~holdflag
  hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [proj] = elproj(pos, method);

% ELPROJ makes a azimuthal projection of a 3D electrode cloud
%  on a plane tangent to the sphere fitted through the electrodes
%  the projection is along the z-axis
%   
%  [proj] = elproj([x, y, z], 'method');
%
% Method should be one of these:
%	  'gnomic'
%	  'stereographic'
%	  'ortographic'
%	  'inverse'
%	  'polar'
%
% Imagine a plane being placed against (tangent to) a globe. If
% a light source inside the globe projects the graticule onto
% the plane the result would be a planar, or azimuthal, map
% projection. If the imaginary light is inside the globe a Gnomonic
% projection results, if the light is antipodal a Sterographic,
% and if at infinity, an Orthographic.
%
% The default projection is a polar projection (BESA like).
% An inverse projection is the opposite of the default polar projection.

% Copyright (C) 2000-2008, Robert Oostenveld
%
% $Log: elproj.m,v $
% Revision 1.4  2008/05/15 10:54:24  roboos
% updated documentation
%
% Revision 1.3  2007/03/20 10:29:35  roboos
% renamed method 'default' into 'polar'
%
% Revision 1.2  2003/03/17 10:37:28  roberto
% improved general help comments and added copyrights
%

x = pos(:,1);
y = pos(:,2);
if size(pos, 2)==3
  z = pos(:,3);
end

if nargin<2
  method='polar';
end

if nargin<3
  secant=1;
end

if strcmp(method, 'orthographic')
  % this method compresses the lowest electrodes very much
  % electrodes on the bottom half of the sphere are folded inwards
  xp = x;
  yp = y;
  num = length(find(z<0));
  str = sprintf('%d electrodes may be folded inwards in orthographic projection\n', num);
  if num
    warning(str);
  end
  proj = [xp yp];

elseif strcmp(method, 'gnomic')
  % the lightsource is in the middle of the sphere
  % electrodes on the equator are projected at infinity
  % electrodes below the equator are not projected at all
  rad = mean(sqrt(x.^2 + y.^2 + z.^2));
  phi = cart2pol(x, y);
  th  = atan(sqrt(x.^2 + y.^2) ./ z);
  xp  = cos(phi) .* tan(th) .* rad;
  yp  = sin(phi) .* tan(th) .* rad;
  num = length(find(th==pi/2 | z<0));
  str = sprintf('removing %d electrodes from gnomic projection\n', num);
  if num
    warning(str);
  end
  xp(find(th==pi/2 | z<0)) = NaN;
  yp(find(th==pi/2 | z<0)) = NaN;
  proj = [xp yp];

elseif strcmp(method, 'stereographic')
  % the lightsource is antipodal (on the south-pole)
  rad = mean(sqrt(x.^2 + y.^2 + z.^2));
  z   = z + rad;
  phi = cart2pol(x, y);
  th  = atan(sqrt(x.^2 + y.^2) ./ z);
  xp  = cos(phi) .* tan(th) .* rad * 2;
  yp  = sin(phi) .* tan(th) .* rad * 2;
  num = length(find(th==pi/2 | z<0));
  str = sprintf('removing %d electrodes from stereographic projection\n', num);
  if num
    warning(str);
  end
  xp(find(th==pi/2 | z<0)) = NaN;
  yp(find(th==pi/2 | z<0)) = NaN;
  proj = [xp yp];
  
elseif strcmp(method, 'inverse')
  % compute the inverse projection of the default angular projection
  [th, r] = cart2pol(x, y);
  [xi, yi, zi] = sph2cart(th, pi/2 - r, 1);
  proj = [xi, yi, zi];

else 
  % use default angular projection
  [az, el, r] = cart2sph(x, y, z);
  [x, y] = pol2cart(az, pi/2 - el);
  proj = [x, y]; 
end

