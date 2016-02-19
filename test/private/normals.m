function [nrm] = normals(pnt, tri, opt)

% NORMALS compute the surface normals of a triangular mesh
% for each triangle or for each vertex
%
% [nrm] = normals(pnt, tri, opt)
% where opt is either 'vertex' or 'triangle'

% Copyright (C) 2002-2007, Robert Oostenveld
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

if nargin<3
  opt='vertex';
elseif (opt(1)=='v' | opt(1)=='V')
  opt='vertex';
elseif (opt(1)=='t' | opt(1)=='T')
  opt='triangle';
else
  error('invalid optional argument');
end

npnt = size(pnt,1);
ntri = size(tri,1);

% shift to center
pnt(:,1) = pnt(:,1)-mean(pnt(:,1),1);
pnt(:,2) = pnt(:,2)-mean(pnt(:,2),1);
pnt(:,3) = pnt(:,3)-mean(pnt(:,3),1);

% compute triangle normals
% nrm_tri = zeros(ntri, 3);
% for i=1:ntri
%   v2 = pnt(tri(i,2),:) - pnt(tri(i,1),:);
%   v3 = pnt(tri(i,3),:) - pnt(tri(i,1),:);
%   nrm_tri(i,:) = cross(v2, v3);
% end

% vectorized version of the previous part
v2 = pnt(tri(:,2),:) - pnt(tri(:,1),:);
v3 = pnt(tri(:,3),:) - pnt(tri(:,1),:);
nrm_tri = cross(v2, v3);


if strcmp(opt, 'vertex')
  % compute vertex normals
  nrm_pnt = zeros(npnt, 3);
  for i=1:ntri
    nrm_pnt(tri(i,1),:) = nrm_pnt(tri(i,1),:) + nrm_tri(i,:);
    nrm_pnt(tri(i,2),:) = nrm_pnt(tri(i,2),:) + nrm_tri(i,:);
    nrm_pnt(tri(i,3),:) = nrm_pnt(tri(i,3),:) + nrm_tri(i,:);
  end
  % normalise the direction vectors to have length one
  nrm = nrm_pnt ./ (sqrt(sum(nrm_pnt.^2, 2)) * ones(1,3));
else
  % normalise the direction vectors to have length one
  nrm = nrm_tri ./ (sqrt(sum(nrm_tri.^2, 2)) * ones(1,3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast cross product to replace the MATLAB standard version
function [c] = cross(a,b)
c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) a(:,1).*b(:,2)-a(:,2).*b(:,1)];

