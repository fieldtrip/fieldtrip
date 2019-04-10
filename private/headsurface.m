function [pos, tri] = headsurface(headmodel, sens, varargin)

% HEADSURFACE constructs a triangulated description of the skin or brain
% surface from a volume conduction model, from a set of electrodes or
% gradiometers, or from a combination of the two. It returns a closed
% surface.
%
% Use as
%   [pos, tri] = headsurface(headmodel, sens, ...)
% where
%   headmodel      = volume conduction model (structure)
%   sens           = electrode or gradiometer array (structure)
%
% Optional arguments should be specified in key-value pairs:
%   surface        = 'skin' or 'brain' (default = 'skin')
%   npos           = number of vertices (default is determined automatic)
%   downwardshift  = boolean, this will shift the lower rim of the helmet down with approximately 1/4th of its radius (default is 1)
%   inwardshift    = number (default = 0)
%   headshape      = string, file containing the head shape

% Copyright (C) 2005-2006, Robert Oostenveld
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

% FIXME perhaps ft_fetch_headshape or ft_prepare_headshape would be a better name for this function

if nargin<1
  headmodel = [];
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
npos          = ft_getopt(varargin, 'npos');                % number of vertices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(headshape)
  % get the surface describing the head shape
  if isstruct(headshape) && isfield(headshape, 'hex')
    headshape = fixpos(headshape);
    fprintf('extracting surface from hexahedral mesh\n');
    headshape = mesh2edge(headshape);
    headshape = poly2tri(headshape);
  elseif isstruct(headshape) && isfield(headshape, 'tet')
    headshape = fixpos(headshape);
    fprintf('extracting surface from tetrahedral mesh\n');
    headshape = mesh2edge(headshape);
  elseif isstruct(headshape) && isfield(headshape, 'tri')
    headshape = fixpos(headshape);
  elseif isstruct(headshape) && isfield(headshape, 'pos')
    headshape = fixpos(headshape);
  elseif isstruct(headshape) && isfield(headshape, 'pnt')
    headshape = fixpos(headshape);
  elseif isnumeric(headshape) && size(headshape,2)==3
    % use the headshape points specified in the configuration
    headshape = struct('pos', headshape);
  elseif ischar(headshape)
    % read the headshape from file
    headshape = ft_read_headshape(headshape);
  end
  if ~isfield(headshape, 'tri')
    for i=1:numel(headshape)
      % generate a closed triangulation from the surface points
      headshape(i).pos = unique(headshape(i).pos, 'rows');
      headshape(i).tri = projecttri(headshape(i).pos);
    end
  end

  % the headshape should be specified as a surface structure with pos and tri
  pos = headshape.pos;
  tri = headshape.tri;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(headmodel) && isfield(headmodel, 'r') && length(headmodel.r)<5
  if length(headmodel.r)==1
    % single sphere model, cannot distinguish between skin and/or brain
    radius = headmodel.r;
    if isfield(headmodel, 'o')
      origin = headmodel.o;
    else
      origin = [0 0 0];
    end
  elseif length(headmodel.r)<5
    % multiple concentric sphere model
    switch surface
      case 'skin'
        % using outermost sphere
        radius = max(headmodel.r);
      case 'brain'
        % using innermost sphere
        radius = min(headmodel.r);
      otherwise
        ft_error('other surfaces cannot be constructed this way');
    end
    if isfield(headmodel, 'o')
      origin = headmodel.o;
    else
      origin = [0 0 0];
    end
  end
  % this requires a specification of the number of vertices
  if isempty(npos)
    npos = 642;
  end
  % construct an evenly tesselated unit sphere
  [pos, tri] = mesh_sphere(npos);

  % scale and translate the vertices
  pos = pos*radius;
  pos(:,1) = pos(:,1) + origin(1);
  pos(:,2) = pos(:,2) + origin(2);
  pos(:,3) = pos(:,3) + origin(3);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ft_headmodeltype(headmodel, 'localspheres')
  % local spheres MEG model, this also requires a gradiometer structure
  grad = sens;
  if ~isfield(grad, 'tra') || ~isfield(grad, 'coilpos')
    ft_error('incorrect specification for the gradiometer array');
  end
  Nchans   = size(grad.tra, 1);
  Ncoils   = size(grad.tra, 2);
  Nspheres = size(headmodel.o, 1);
  if Nspheres~=Ncoils
    ft_error('there should be just as many spheres as coils');
  end
  % for each coil, determine a surface point using the corresponding sphere
  vec = grad.coilpos - headmodel.o;
  nrm = sqrt(sum(vec.^2,2));
  vec = vec ./ [nrm nrm nrm];
  pos = headmodel.o + vec .* [headmodel.r headmodel.r headmodel.r];
  pos = unique(pos, 'rows');
  %  make a triangularization that is needed to find the rim of the helmet
  prj = elproj(pos);
  tri = delaunay(prj(:,1),prj(:,2));
  % find the lower rim of the helmet shape
  [pos, line] = find_mesh_edge(pos, tri);
  edgeind     = unique(line(:));
  % shift the lower rim of the helmet shape down with approximately 1/4th of its radius
  if downwardshift
    % determine the extent of the volume conduction model
    dist = mean(sqrt(sum((pos - repmat(mean(pos,1), size(pos,1), 1)).^2, 2)));
    dist = dist/4;
    pos(edgeind,3) = pos(edgeind,3) - dist;
  end
  % construct the triangulation of the surface
  tri = projecttri(pos);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ft_headmodeltype(headmodel, 'bem') ||  ft_headmodeltype(headmodel, 'singleshell')
  % volume conduction model with triangulated boundaries
  switch surface
    case 'skin'
      if ~isfield(headmodel, 'skin')
        headmodel.skin   = find_outermost_boundary(headmodel.bnd);
      end
      pos = headmodel.bnd(headmodel.skin).pos;
      tri = headmodel.bnd(headmodel.skin).tri;
    case 'brain'
      if ~isfield(headmodel, 'source')
        headmodel.source  = find_innermost_boundary(headmodel.bnd);
      end
      pos = headmodel.bnd(headmodel.source).pos;
      tri = headmodel.bnd(headmodel.source).tri;
    otherwise
      ft_error('other surfaces cannot be constructed this way');
  end
end

% retriangulate the skin/brain/cortex surface to the desired number of vertices
if ~isempty(npos) && size(pos,1)~=npos
  [pnt2, tri2] = mesh_sphere(npos);
  [pos, tri]   = retriangulate(pos, tri, pnt2, tri2, 2);
end

% shift the surface inward with a certain amount
if ~isempty(inwardshift) && inwardshift~=0
  ori = normals(pos, tri, 'vertex');
  % FIXME in case of a icosahedron projected onto a localspheres model, the
  % surfaceorientation for the lower rim points fails, causing problems
  % with the inward shift
  tmp = surfaceorientation(pos, tri, ori);
  % the orientation of the normals should be pointing to the outside of the surface
  if tmp==1
    % the normals are outward oriented
    % nothing to do
  elseif tmp==-1
    % the normals are inward oriented
    tri = fliplr(tri);
    ori = -ori;
  else
    ft_warning('cannot determine the orientation of the vertex normals');
    % nothing to do
  end
  % the orientation is outward, hence shift with a negative amount
  pos = pos - inwardshift * ori;
end
