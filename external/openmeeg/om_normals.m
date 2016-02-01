function [nrm] = normals(pos, tri, opt)

% NORMALS compute the surface normals of a triangular mesh
% for each triangle or for each vertex
%
% Use as
%   nrm = normals(pos, tri, opt)
% where opt is either 'vertex' or 'triangle'

% Copyright (C) 2002-2007, Robert Oostenveld

if nargin<3
  opt='vertex';
elseif (opt(1)=='v' | opt(1)=='V')
  opt='vertex';
elseif (opt(1)=='t' | opt(1)=='T')
  opt='triangle';
else
  error('invalid optional argument');
end

npos = size(pos,1);
ntri = size(tri,1);

% shift to center
pos(:,1) = pos(:,1)-mean(pos(:,1),1);
pos(:,2) = pos(:,2)-mean(pos(:,2),1);
pos(:,3) = pos(:,3)-mean(pos(:,3),1);

% compute triangle normals
nrm_tri = zeros(ntri, 3);
for i=1:ntri
  v2 = pos(tri(i,2),:) - pos(tri(i,1),:);
  v3 = pos(tri(i,3),:) - pos(tri(i,1),:);
  nrm_tri(i,:) = cross(v2, v3);
end

if strcmp(opt, 'vertex')
  % compute vertex normals
  nrm_pos = zeros(npos, 3);
  for i=1:ntri
    nrm_pos(tri(i,1),:) = nrm_pos(tri(i,1),:) + nrm_tri(i,:);
    nrm_pos(tri(i,2),:) = nrm_pos(tri(i,2),:) + nrm_tri(i,:);
    nrm_pos(tri(i,3),:) = nrm_pos(tri(i,3),:) + nrm_tri(i,:);
  end
  % normalise the direction vectors to have length one
  nrm = nrm_pos ./ (sqrt(sum(nrm_pos.^2, 2)) * ones(1,3));
else
  % normalise the direction vectors to have length one
  nrm = nrm_tri ./ (sqrt(sum(nrm_tri.^2, 2)) * ones(1,3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast cross product to replace the Matlab standard version
function [c] = cross(a,b)
c = [a(2)*b(3)-a(3)*b(2) a(3)*b(1)-a(1)*b(3) a(1)*b(2)-a(2)*b(1)];

