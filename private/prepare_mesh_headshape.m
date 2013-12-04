function bnd = prepare_mesh_headshape(cfg)

% PREPARE_MESH_HEADSHAPE
%
% See also PREPARE_MESH_MANUAL, PREPARE_MESH_SEGMENTATION

% Copyrights (C) 2009, Robert Oostenveld
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

% get the specific options
cfg.headshape    = ft_getopt(cfg, 'headshape');

if isa(cfg.headshape, 'config')
  % convert the nested config-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

% get the surface describing the head shape
if isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
  % use the headshape surface specified in the configuration
  headshape = cfg.headshape;
elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
  % use the headshape points specified in the configuration
  headshape.pnt = cfg.headshape;
elseif ischar(cfg.headshape)
  % read the headshape from file
  headshape = ft_read_headshape(cfg.headshape);
else
  error('cfg.headshape is not specified correctly')
end

% usually a headshape only describes a single surface boundaries, but there are cases
% that multiple surfaces are included, e.g. skin_surface, outer_skull_surface, inner_skull_surface
nbnd = numel(headshape);

if ~isfield(headshape, 'tri')
  % generate a closed triangulation from the surface points
  for i=1:nbnd
    headshape(i).pnt = unique(headshape(i).pnt, 'rows');
    headshape(i).tri = projecttri(headshape(i).pnt);
  end
end

if ~isempty(cfg.numvertices) && ~strcmp(cfg.numvertices, 'same')
  for i=1:nbnd
    tri1 = headshape(i).tri;
    pnt1 = headshape(i).pnt;
    % The number of vertices is multiplied by 3 in order to have more
    % points on the original mesh than on the sphere mesh (see below).
    % The rationale for this is that every projection point on the sphere
    % has three corresponding points on the mesh
    if (cfg.numvertices>size(pnt1,1))
      [tri1, pnt1] = refinepatch(headshape(i).tri, headshape(i).pnt, 3*cfg.numvertices);
    else
      [tri1, pnt1] = reducepatch(headshape(i).tri, headshape(i).pnt, 3*cfg.numvertices);
    end
    
    % remove double vertices
    [pnt1, tri1] = remove_double_vertices(pnt1, tri1);
    
    % replace the probably unevenly distributed triangulation with a regular one
    % and retriangulate it to the desired accuracy
    [pnt2, tri2] = mysphere(cfg.numvertices); % this is a regular triangulation
    [pnt1, tri1] = retriangulate(pnt1, tri1, pnt2, tri2, 2);
    [pnt1, tri1] = fairsurface(pnt1, tri1, 1);% this helps redistribute the superimposed points
    
    % remove double vertices
    [headshape(i).pnt,headshape(i).tri] = remove_double_vertices(pnt1, tri1);
    fprintf('returning %d vertices, %d triangles\n', size(headshape(i).pnt,1), size(headshape(i).tri,1));
  end
end

% the output should only describe one or multiple boundaries and should not
% include any other fields
bnd = rmfield(headshape, setdiff(fieldnames(headshape), {'pnt', 'tri'}));

function [tri1, pnt1] = refinepatch(tri, pnt, numvertices)
fprintf('the original mesh has %d vertices against the %d requested\n',size(pnt,1),numvertices/3);
fprintf('trying to refine the compartment...\n');
[pnt1, tri1] = refine(pnt, tri, 'updown', numvertices);

function [pnt, tri] = mysphere(N)
% This is a copy of MSPHERE without the confusing output message
% Returns a triangulated sphere with approximately M vertices
% that are nicely distributed over the sphere. The vertices are aligned
% along equally spaced horizontal contours according to an algorithm of
% Dave Russel.
%
% Use as
%  [pnt, tri] = msphere(M)
%
% See also SPHERE, NSPHERE, ICOSAHEDRON, REFINE
% Copyright (C) 1994, Dave Rusin

storeM    = [];
storelen  = [];
increaseM = 0;
while (1)
  
  % put a single vertex at the top
  phi = [0];
  th  = [0];
  
  M = round((pi/4)*sqrt(N)) + increaseM;
  for k=1:M
    newphi = (k/M)*pi;
    Q = round(2*M*sin(newphi));
    for j=1:Q
      phi(end+1) = newphi;
      th(end+1)  = (j/Q)*2*pi;
      % in case of even number of contours
      if mod(M,2) & k>(M/2)
        th(end) = th(end) + pi/Q;
      end
    end
  end
  
  % put a single vertex at the bottom
  phi(end+1) = [pi];
  th(end+1)  = [0];
  
  % store this vertex packing
  storeM(end+1).th  = th;
  storeM(end  ).phi = phi;
  storelen(end+1) = length(phi);
  if storelen(end)>N
    break;
  else
    increaseM = increaseM+1;
    % fprintf('increasing M by %d\n', increaseM);
  end
end

% take the vertex packing that most closely matches the requirement
[m, i] = min(abs(storelen-N));
th  = storeM(i).th;
phi = storeM(i).phi;

% convert from spherical to cartehsian coordinates
[x, y, z] = sph2cart(th, pi/2-phi, 1);
pnt = [x' y' z'];
tri = convhulln(pnt);

function [pntR, triR] = remove_double_vertices(pnt, tri)

% REMOVE_VERTICES removes specified vertices from a triangular mesh
% renumbering the vertex-indices for the triangles and removing all
% triangles with one of the specified vertices.
%
% Use as
%   [pnt, tri] = remove_double_vertices(pnt, tri)

pnt1 = unique(pnt, 'rows');
keeppnt   = find(ismember(pnt1,pnt,'rows'));
removepnt = setdiff([1:size(pnt,1)],keeppnt);

npnt = size(pnt,1);
ntri = size(tri,1);

if all(removepnt==0 | removepnt==1)
  removepnt = find(removepnt);
end

% remove the vertices and determine the new numbering (indices) in numb
keeppnt = setdiff(1:npnt, removepnt);
numb    = zeros(1,npnt);
numb(keeppnt) = 1:length(keeppnt);

% look for triangles referring to removed vertices
removetri = false(ntri,1);
removetri(ismember(tri(:,1), removepnt)) = true;
removetri(ismember(tri(:,2), removepnt)) = true;
removetri(ismember(tri(:,3), removepnt)) = true;

% remove the vertices and triangles
pntR = pnt(keeppnt, :);
triR = tri(~removetri,:);

% renumber the vertex indices for the triangles
triR = numb(triR);

function [pnt1, tri1] = fairsurface(pnt, tri, N)

% FAIRSURFACE modify the mesh in order to reduce overlong edges, and
% smooth out "rough" areas. This is a non-shrinking smoothing algorithm.
% The procedure uses an elastic model : At each vertex, the neighbouring
% triangles and vertices connected directly are used. Each edge is
% considered elastic and can be lengthened or shortened, depending
% on their length. Displacement are done in 3D, so that holes and
% bumps are attenuated.
%
% Use as
%   [pnt, tri] = fairsurface(pnt, tri, N);
% where N is the number of smoothing iterations.
%
% This implements:
%   G.Taubin, A signal processing approach to fair surface design, 1995

% This function corresponds to spm_eeg_inv_ElastM
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%                    Christophe Phillips & Jeremie Mattout
% spm_eeg_inv_ElastM.m 1437 2008-04-17 10:34:39Z christophe
%
% $Id$

ts = [];
ts.XYZmm = pnt';
ts.tri   = tri';
ts.nr(1) = size(pnt,1);
ts.nr(2) = size(tri,1);

% Connection vertex-to-vertex
%--------------------------------------------------------------------------
M_con = sparse([ts.tri(1,:)';ts.tri(1,:)';ts.tri(2,:)';ts.tri(3,:)';ts.tri(2,:)';ts.tri(3,:)'], ...
  [ts.tri(2,:)';ts.tri(3,:)';ts.tri(1,:)';ts.tri(1,:)';ts.tri(3,:)';ts.tri(2,:)'], ...
  ones(ts.nr(2)*6,1),ts.nr(1),ts.nr(1));

kpb   = .1;                       % Cutt-off frequency
lam   = .5; mu = lam/(lam*kpb-1); % Parameters for elasticity.
XYZmm = ts.XYZmm;

% smoothing iterations
%--------------------------------------------------------------------------
for j=1:N
  
  XYZmm_o = zeros(3,ts.nr(1)) ;
  XYZmm_o2 = zeros(3,ts.nr(1)) ;
  
  for i=1:ts.nr(1)
    ln = find(M_con(:,i));
    d_i = sqrt(sum((XYZmm(:,ln)-XYZmm(:,i)*ones(1,length(ln))).^2));
    if sum(d_i)==0
      w_i = zeros(size(d_i));
    else
      w_i = d_i/sum(d_i);
    end
    XYZmm_o(:,i) = XYZmm(:,i) + ...
      lam * sum((XYZmm(:,ln)-XYZmm(:,i)*ones(1,length(ln))).*(ones(3,1)*w_i),2);
  end
  
  for i=1:ts.nr(1)
    ln = find(M_con(:,i));
    d_i = sqrt(sum((XYZmm(:,ln)-XYZmm(:,i)*ones(1,length(ln))).^2));
    if sum(d_i)==0
      w_i = zeros(size(d_i));
    else
      w_i = d_i/sum(d_i);
    end
    XYZmm_o2(:,i) = XYZmm_o(:,i) + ...
      mu * sum((XYZmm_o(:,ln)-XYZmm_o(:,i)*ones(1,length(ln))).*(ones(3,1)*w_i),2);
  end
  
  XYZmm = XYZmm_o2;
  
end

% collect output results
%--------------------------------------------------------------------------

pnt1 = XYZmm';
tri1 = tri;

if 0
  % this is some test/demo code
  bnd = [];
  [bnd.pnt, bnd.tri] = icosahedron162;
  
  scale = 1+0.3*randn(size(pnt,1),1);
  bnd.pnt = bnd.pnt .* [scale scale scale];
  
  figure
  ft_plot_mesh(bnd)
  
  [bnd.pnt, bnd.tri] = fairsurface(bnd.pnt, bnd.tri, 10);
  
  figure
  ft_plot_mesh(bnd)
end
