function [depth] = sourcedepth(pos, vol)

% SOURCEDEPTH computes the distance from the source to the surface of
% the source compartment (usually the brain).
%
% Use as
%   depth = sourcedepth(pos, vol);
% where
%   pos     Nx3 matrix with the position of N sources
%   vol     structure describing volume condition model
%
% A negative depth indicates that the source is inside the source
% compartment, positive indicates outside.
%
% See also FIND_INSIDE_VOL

% Copyright (C) 2007-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% determine the type of volume conduction model
switch voltype(vol)

% single-sphere or multiple concentric spheres
case {'singlesphere', 'concentric'}
  if ~isfield(vol, 'source')
    % locate the innermost compartment and remember it
    [dum, vol.source] = min(vol.r);
  end
  if isfield(vol, 'o')
    % shift dipole positions toward origin of sphere
    tmp = pos - repmat(vol.o, size(pos,1), 1);
  else
    tmp = pos;
  end
  depth = sqrt(sum(tmp.^2, 2))-vol.r(vol.source); % positive if outside, negative if inside

% boundary element model
case {'bem' 'dipoli', 'bemcp', 'asa', 'avo', 'nolte', 'neuromag'}
  if isfield(vol, 'source')
    % use the specified source compartment
    pnt = vol.bnd(vol.source).pnt;
    tri = vol.bnd(vol.source).tri;
  else
    % locate the innermost compartment and remember it
    vol.source = find_innermost_boundary(vol.bnd);
    pnt = vol.bnd(vol.source).pnt;
    tri = vol.bnd(vol.source).tri;
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

