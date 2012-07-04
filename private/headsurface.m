function [pnt, tri] = headsurface(vol, sens, varargin);

% HEADSURFACE constructs a triangulated description of the skin or brain
% surface from a volume conduction model, from a set of electrodes or
% gradiometers, or from a combination of the two. It returns a closed
% surface.
%
% Use as
%   [pnt, tri] = headsurface(vol, sens, ...)
% where
%   vol            volume conduction model (structure)
%   sens           electrode or gradiometer array (structure)
%
% Optional arguments should be specified in key-value pairs:
%   surface        = 'skin' or 'brain' (default = 'skin')
%   npnt           = number of vertices (default is determined automatic)
%   downwardshift  = boolean, this will shift the lower rim of the helmet down with approximately 1/4th of its radius (default is 1)
%   inwardshift    = number (default = 0)
%   headshape      = string, file containing the head shape

% Copyright (C) 2005-2006, Robert Oostenveld
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

if nargin<1
  vol = [];
end

if nargin<2
  sens = [];
end

if ~isempty(sens)
  sens = ft_datatype_sens(sens);
end

if nargin<3
  varargin = {};
end

% parse the optional input arguments
surface       = ft_getopt(varargin, 'surface', 'skin');     % skin or brain
downwardshift = ft_getopt(varargin, 'downwardshift', true); % boolean
inwardshift   = ft_getopt(varargin, 'inwardshift');         % number
headshape     = ft_getopt(varargin, 'headshape');           % CTF *.shape file
npnt          = ft_getopt(varargin, 'npnt');                % number of vertices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(headshape)
  % the headshape should be specified as a surface structure with pnt and tri
  pnt = headshape.pnt;
  tri = headshape.tri;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(vol) && isfield(vol, 'r') && length(vol.r)<5
  if length(vol.r)==1
    % single sphere model, cannot distinguish between skin and/or brain
    radius = vol.r;
    if isfield(vol, 'o')
      origin = vol.o;
    else
      origin = [0 0 0];
    end
  elseif length(vol.r)<5
    % multiple concentric sphere model
    switch surface
      case 'skin'
        % using outermost sphere
        radius = max(vol.r);
      case 'brain'
        % using innermost sphere
        radius = min(vol.r);
      otherwise
        error('other surfaces cannot be constructed this way');
    end
    if isfield(vol, 'o')
      origin = vol.o;
    else
      origin = [0 0 0];
    end
  end
  % this requires a specification of the number of vertices
  if isempty(npnt)
    npnt = 642;
  end
  % construct an evenly tesselated unit sphere
  switch npnt
    case 2562
      [pnt, tri] = icosahedron2562;
    case 642
      [pnt, tri] = icosahedron642;
    case 162
      [pnt, tri] = icosahedron162;
    case 42
      [pnt, tri] = icosahedron42;
    case 12
      [pnt, tri] = icosahedron;
    otherwise
      [pnt, tri] = ksphere(npnt);
  end
  % scale and translate the vertices
  pnt = pnt*radius;
  pnt(:,1) = pnt(:,1) + origin(1);
  pnt(:,2) = pnt(:,2) + origin(2);
  pnt(:,3) = pnt(:,3) + origin(3);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ft_voltype(vol, 'localspheres')
  % local spheres MEG model, this also requires a gradiometer structure
  grad = sens;
  if ~isfield(grad, 'tra') || ~isfield(grad, 'coilpos')
    error('incorrect specification for the gradiometer array');
  end
  Nchans   = size(grad.tra, 1);
  Ncoils   = size(grad.tra, 2);
  Nspheres = size(vol.o, 1);
  if Nspheres~=Ncoils
    error('there should be just as many spheres as coils');
  end
  % for each coil, determine a surface point using the corresponding sphere
  vec = grad.coilpos - vol.o;
  nrm = sqrt(sum(vec.^2,2));
  vec = vec ./ [nrm nrm nrm];
  pnt = vol.o + vec .* [vol.r vol.r vol.r];
  pnt = unique(pnt, 'rows');
  %  make a triangularization that is needed to find the rim of the helmet
  prj = elproj(pnt);
  tri = delaunay(prj(:,1),prj(:,2));
  % find the lower rim of the helmet shape
  [pnt, line] = find_mesh_edge(pnt, tri);
  edgeind     = unique(line(:));
  % shift the lower rim of the helmet shape down with approximately 1/4th of its radius
  if downwardshift
    % determine the extent of the volume conduction model
    dist = mean(sqrt(sum((pnt - repmat(mean(pnt,1), size(pnt,1), 1)).^2, 2)));
    dist = dist/4;
    pnt(edgeind,3) = pnt(edgeind,3) - dist;
  end
  % construct the triangulation of the surface
  tri = projecttri(pnt);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ft_voltype(vol, 'bem') ||  ft_voltype(vol, 'nolte')
  % volume conduction model with triangulated boundaries
  switch surface
    case 'skin'
      if ~isfield(vol, 'skin')
        vol.skin   = find_outermost_boundary(vol.bnd);
      end
      pnt = vol.bnd(vol.skin).pnt;
      tri = vol.bnd(vol.skin).tri;
    case 'brain'
      if ~isfield(vol, 'source')
        vol.source  = find_innermost_boundary(vol.bnd);
      end
      pnt = vol.bnd(vol.source).pnt;
      tri = vol.bnd(vol.source).tri;
    otherwise
      error('other surfaces cannot be constructed this way');
  end
end

% retriangulate the skin/brain/cortex surface to the desired number of vertices
if ~isempty(npnt) && size(pnt,1)~=npnt
  switch npnt
    case 2562
      [pnt2, tri2] = icosahedron2562;
    case 642
      [pnt2, tri2] = icosahedron642;
    case 162
      [pnt2, tri2] = icosahedron162;
    case 42
      [pnt2, tri2] = icosahedron42;
    case 12
      [pnt2, tri2] = icosahedron;
    otherwise
      [pnt2, tri2] = ksphere(npnt);
  end
  [pnt, tri] = retriangulate(pnt, tri, pnt2, tri2, 2);
end

% shift the surface inward with a certain amount
if ~isempty(inwardshift) && inwardshift~=0
  ori = normals(pnt, tri, 'vertex');
  % FIXME in case of a icosahedron projected onto a localspheres model, the
  % surfaceorientation for th elower rim points fails, causing problems
  % with the inward shift
  tmp = surfaceorientation(pnt, tri, ori);
  % the orientation of the normals should be pointing to the outside of the surface
  if tmp==1
    % the normals are outward oriented
    % nothing to do
  elseif tmp==-1
    % the normals are inward oriented
    warning('the normals of the surface triangulation are inward oriented');
    tri = fliplr(tri);
    ori = -ori;
  else
    warning('cannot determine the orientation of the vertex normals');
    % nothing to do
  end
  % the orientation is outward, hence shift with a negative amount
  pnt = pnt - inwardshift * ori;
end
