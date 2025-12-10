function nv = nodesurfnorm(node, elem)
%  nv=nodesurfnorm(node,elem)
%
%  calculate a nodal norm for each vertix on a surface mesh (surface
%   can only be triangular or cubic)
%
%  author: Qianqian Fang <q.fang at neu.edu>
%  date: 12/12/2008
%
% parameters:
%      node: node coordinate of the surface mesh (nn x 3)
%      elem: element list of the surface mesh (3 columns for
%            triangular mesh, 4 columns for cubic surface mesh)
%      pt: points to be projected, 3 columns for x,y and z respectively
%
% outputs:
%      nv: nodal norms (vector) calculated from nodesurfnorm.m
%          with dimensions of (size(v,1),3)
%
% Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%
% this function is part of "metch" toobox, see COPYING for license

nn = size(node, 1);
ne = size(elem, 1);
nedim = size(elem, 2);

ev = surfacenorm(node, elem);

nv = zeros(nn, 3);
ev2 = repmat(ev, 1, 3);
for i = 1:ne
    nv(elem(i, :), :) = nv(elem(i, :), :) + reshape(ev2(i, :), 3, 3)';
end
nvnorm = sqrt(sum(nv .* nv, 2));
idx = find(nvnorm > 0);
if (length(idx) < nn)
    warning(['found interior nodes, their norms will be set to zeros; to remove ', ...
             'them, please use removeisolatednodes.m from iso2mesh toolbox']);

    nv(idx, :) = nv(idx, :) ./ repmat(nvnorm(idx), 1, 3);
else
    nv = nv ./ repmat(nvnorm, 1, 3);
end
