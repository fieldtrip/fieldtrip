function [depth] = ft_sourcedepth(dippos, headmodel)

% FT_SOURCEDEPTH computes the distance from the source to the surface of
% the source compartment (usually the brain) in the volume conduction model.
%
% Use as
%   depth = ft_sourcedepth(dippos, headmodel);
% where
%   dippos    =  Nx3 matrix with the position of N sources
%   headmodel =  structure describing volume condition model
%
% A negative depth indicates that the source is inside the source
% compartment, positive indicates outside.
%
% See also FIND_INSIDE_VOL

% Copyright (C) 2007-2008, Robert Oostenveld
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

% ensure the representation of the headmodel to be up-to-date
headmodel = ft_datatype_headmodel(headmodel);

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
    tmp = dippos - repmat(headmodel.o, size(dippos,1), 1);
  else
    tmp = dippos;
  end
  depth = sqrt(sum(tmp.^2, 2))-headmodel.r(headmodel.source); % positive if outside, negative if inside

% boundary element model
case {'bem' 'dipoli', 'bemcp', 'asa', 'singleshell', 'neuromag','openmeeg'}
  if isfield(headmodel, 'source')
    % use the specified source compartment
    pos = headmodel.bnd(headmodel.source).pos;
    tri = headmodel.bnd(headmodel.source).tri;
  else
    % locate the innermost compartment and remember it
    headmodel.source = find_innermost_boundary(headmodel.bnd);
    pos = headmodel.bnd(headmodel.source).pos;
    tri = headmodel.bnd(headmodel.source).tri;
  end
  inside = bounding_mesh(dippos, pos, tri);
  ntri   = size(tri,1);
  npos   = size(dippos,1);
  dist   = zeros(ntri, 1);
  depth  = zeros(npos, 1);
  for i=1:npos
    for j=1:ntri
      v1 = pos(tri(j,1),:);
      v2 = pos(tri(j,2),:);
      v3 = pos(tri(j,3),:);
      [proj, dist(j)] = ptriproj(v1, v2, v3, dippos(i,:), 1);
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

