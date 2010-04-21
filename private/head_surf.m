function [pnt,tri]=head_surf(vol,grad,flag);

% HEAD_SURF determines the head surface from a multisphere headmodel
% and a set of coils in the MEG gradiometer system and returns a
% triangulation
% 
% Use as
%   [pnt, tri] = head_surf(vol, grad, flag) 
% where
%   grad    gradiometer definition
%   vol     multisphere volume conductor definition
%
% If flag=1 the lower rim of the helmet-shaped head surface will
% be shifted downward, if flag=0 it will not be shifted downward.
% The default is flag=1.
%
% The head surface triangulation can be used in combination with
% find_inside_vol to determine which voxels of a beamformer scan 
% are located within the head (region of interest).

% Copyright (C) Jan-Matthijs Schoffelen
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

if nargin<3
  flag = 1;
end

Nchans = size(grad.tra, 1);
Ncoils = size(grad.tra, 2);

% for each coil, determine a surface point using the corresponding sphere
vec = grad.pnt - vol.o;
nrm = sqrt(sum(vec.^2,2));
vec = vec ./ [nrm nrm nrm];
pnt = vol.o + vec .* [vol.r vol.r vol.r];

%  make a triangularization that is needed to find the rim of the helmet
prj = elproj(pnt);
tri = delaunay(prj(:,1),prj(:,2));

% find the lower rim of the helmet shape
[pnt,line] = find_mesh_edge(pnt, tri);
edgeind    = unique(line(:));

if flag
  % shift the lower rim of the helmet shape down with approximately 1/4th of its radius
  dist = mean(sqrt(sum((pnt - repmat(mean(pnt,1), Ncoils, 1)).^2, 2)));
  dist = dist/4;
  pnt(edgeind,3) = pnt(edgeind,3) - dist; 
end

% use matlab triangulation algorithm to determine convex hull, this is the
% final triangulation which makes a nice headshape
tri = convhulln(pnt);
