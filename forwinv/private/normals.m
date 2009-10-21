function [nrm] = normals(pnt, dhk, opt);

% NORMALS compute the surface normals of a triangular mesh
% for each triangle or for each vertex
%
% [nrm] = normals(pnt, dhk, opt)
% where opt is either 'vertex' or 'triangle'

% Copyright (C) 2002-2007, Robert Oostenveld
%
% $Log: normals.m,v $
% Revision 1.5  2008/10/21 09:00:34  roboos
% use local version of "cross" function to speed it up
%
% Revision 1.4  2007/07/26 07:13:08  roboos
% shift center of vertices to origin
%
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

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
ndhk = size(dhk,1);

% shift to center
pnt(:,1) = pnt(:,1)-mean(pnt(:,1),1);
pnt(:,2) = pnt(:,2)-mean(pnt(:,2),1);
pnt(:,3) = pnt(:,3)-mean(pnt(:,3),1);

% compute triangle normals
nrm_dhk = zeros(ndhk, 3);
for i=1:ndhk
  v2 = pnt(dhk(i,2),:) - pnt(dhk(i,1),:);
  v3 = pnt(dhk(i,3),:) - pnt(dhk(i,1),:);
  nrm_dhk(i,:) = cross(v2, v3);
end

if strcmp(opt, 'vertex')
  % compute vertex normals
  nrm_pnt = zeros(npnt, 3);
  for i=1:ndhk
    nrm_pnt(dhk(i,1),:) = nrm_pnt(dhk(i,1),:) + nrm_dhk(i,:);
    nrm_pnt(dhk(i,2),:) = nrm_pnt(dhk(i,2),:) + nrm_dhk(i,:);
    nrm_pnt(dhk(i,3),:) = nrm_pnt(dhk(i,3),:) + nrm_dhk(i,:);
  end
  % normalise the direction vectors to have length one
  nrm = nrm_pnt ./ (sqrt(sum(nrm_pnt.^2, 2)) * ones(1,3));
else
  % normalise the direction vectors to have length one
  nrm = nrm_dhk ./ (sqrt(sum(nrm_dhk.^2, 2)) * ones(1,3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast cross product to replace the Matlab standard version
function [c] = cross(a,b)
c = [a(2)*b(3)-a(3)*b(2) a(3)*b(1)-a(1)*b(3) a(1)*b(2)-a(2)*b(1)];

