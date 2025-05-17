function [d2surf, cn] = dist2surf(node, nv, p, cn)
%  [d2surf,cn]=dist2surf(node,nv,p)
%
%  calculate the distances from a point cloud to a surface, and return
%  the indices of the closest surface node
%
%  author: Qianqian Fang <q.fang at neu.edu>
%  date: 12/12/2008
%
% parameters:
%      node: node coordinate of the surface mesh (nn x 3)
%      nv: nodal norms (vector) calculated from nodesurfnorm.m
%          with dimensions of (size(node,1),3), this can be
%          calcuated from nodesurfnorm.m
%      pt: points to be calculated, 3 columns for x,y and z respectively
%
% outputs:
%      d2surf: a vector of length of p, the distances from p(i) to the surface
%      cn: a integer vector with the length of p, the indices of the closest surface node
%
% Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%
% this function is part of "metch" toobox, see COPYING for license

if (nargin < 4)
    nn = size(node, 1);
    pnum = size(p, 1);
    mindist = zeros(pnum, 1);
    cn = zeros(pnum, 1);
    for i = 1:pnum
        d0 = node - repmat(p(i, :), nn, 1);
        d0 = sum(d0 .* d0, 2);
        [mindist(i), cn(i)] = min(d0);
    end
end
d2surf = abs(sum(nv(cn, :) .* (p - node(cn, :)), 2));
