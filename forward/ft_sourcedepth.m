function [depth] = ft_sourcedepth(pos, headmodel)

% FT_SOURCEDEPTH computes the distance from the source to the surface of
% the source compartment (usually the brain) in the volume conduction model.
%
% Use as
%   depth = ft_sourcedepth(pos, headmodel);
% where
%   pos       =  Nx3 matrix with the position of N sources
%   headmodel =  structure describing volume condition model
%
% A negative depth indicates that the source is inside the source
% compartment, positive indicates outside.
%
% See also FIND_INSIDE_VOL

% Copyright (C) 2007-2008, Robert Oostenveld
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

% determine the type of volume conduction model
switch ft_voltype(headmodel)

% single-sphere or multiple concentric spheres
case {'singlesphere', 'concentricspheres'}
  if ~isfield(headmodel, 'source')
    % locate the innermost compartment and remember it
    [dum, headmodel.source] = min(headmodel.r);
  end
  if isfield(headmodel, 'o')
    % shift dipole positions toward origin of sphere
    tmp = pos - repmat(headmodel.o, size(pos,1), 1);
  else
    tmp = pos;
  end
  depth = sqrt(sum(tmp.^2, 2))-headmodel.r(headmodel.source); % positive if outside, negative if inside

% boundary element model
case {'bem' 'dipoli', 'bemcp', 'asa', 'singleshell', 'neuromag','openmeeg'}
  if isfield(headmodel, 'source')
    % use the specified source compartment
    pnt = headmodel.bnd(headmodel.source).pnt;
    tri = headmodel.bnd(headmodel.source).tri;
  else
    % locate the innermost compartment and remember it
    headmodel.source = find_innermost_boundary(headmodel.bnd);
    pnt = headmodel.bnd(headmodel.source).pnt;
    tri = headmodel.bnd(headmodel.source).tri;
  end
  inside = bounding_mesh(pos, pnt, tri);
  ntri   = size(tri,1);
  npos   = size(pos,1);
  dist   = zeros(ntri, 1);
  depth  = zeros(npos, 1);
  for i=1:npos
    for j=1:ntri
      v1 = pnt(tri(j,1),:);
      v2 = pnt(tri(j,2),:);
      v3 = pnt(tri(j,3),:);
      [proj, dist(j)] = ptriproj(v1, v2, v3, pos(i,:), 1);
    end
    if inside(i)
      depth(i) = -min(dist);
    else
      depth(i) = min(dist);
    end
  end

% unsupported volume conductor model
otherwise
  error('upsupported volume conductor model');
end

