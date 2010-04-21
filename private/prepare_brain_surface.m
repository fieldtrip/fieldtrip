function [pnt] = prepare_brain_surface(cfg, grad, vol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that makes a brain surface from the headhape and places a lot
% of dipoles on it. The dipoles subsequently can be used for a simplified 
% distributed source model.
%
% This function uses the following fields from the configuration
%   cfg.headshape
%   cfg.spheremesh
%   cfg.inwardshift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2004, Robet Oostenveld
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

% create a head surface, which by default is assumed to correspond to the skin
if strcmp(cfg.headshape, 'headmodel')
  if length(vol.r)==1
    % single sphere model, use a refined icosahedron 
    fprintf('constructing brain surface from single sphere model\n');
    if cfg.spheremesh==2562
      [pnt, tri] = icosahedron2562;
    elseif cfg.spheremesh==642
      [pnt, tri] = icosahedron642;
    elseif cfg.spheremesh==162
      [pnt, tri] = icosahedron162;
    elseif cfg.spheremesh==42
      [pnt, tri] = icosahedron42;
    elseif cfg.spheremesh==12
      [pnt, tri] = icosahedron;
    else
      [pnt, tri] = ksphere(cfg.spheremesh);
    end
    % scale the sourcemodel sphere to the size of the headmodel sphere
    pnt = pnt * vol.r;
    pnt(:,1) = pnt(:,1) + vol.o(1);
    pnt(:,2) = pnt(:,2) + vol.o(2);
    pnt(:,3) = pnt(:,3) + vol.o(3);
  else
    % multiple sphere model, use the points on the skin surface
    fprintf('constructing brain surface from multiple sphere model\n');
    [pnt, tri] = head_surf(vol, grad, 0);
    prj = elproj(pnt);
    tri = delaunay(prj(:,1),prj(:,2));
  end
else
  fprintf('constructing brain surface from headshape file\n');
  % read the headshape from file
  shape = read_headshape(cfg.headshape);
  pnt = shape.pnt;
  prj = elproj(pnt);
  tri = delaunay(prj(:,1), prj(:,2));
  % the number of triangles is approximately twice the number of vertices
  [tri, pnt] = reducepatch(tri, pnt, 2*cfg.spheremesh);
end

% correct the orientation of the triangles
sel=find(solid_angle(pnt,tri)<0);
tri(sel,:) = fliplr(tri(sel,:));

% shift the head surface inward with a certain amount
ori = normals(pnt, tri, 'vertex');
pnt = pnt + cfg.inwardshift * ori;
