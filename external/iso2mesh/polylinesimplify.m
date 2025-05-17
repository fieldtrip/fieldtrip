function [newnodes, len] = polylinesimplify(nodes, minangle)
%
% [newnodes, len]=polylinesimplify(nodes, minangle)
%
% Calculate a simplified polyline by removing nodes where two adjacent
% segment have an angle less than a specified limit
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    node: an N x 3 array defining each vertex of the polyline in
%          sequential order
%    minangle:(optional) minimum segment angle in radian, if not given, use
%          0.75*pi
%
% output:
%    newnodes: the updated node list; start/end will not be removed
%    len: the length of each segment between the start and the end points
%
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v3 or later, see LICENSE.txt for details
%

if (nargin < 2)
    minangle = 0.75 * pi;
end

v = segvec(nodes(1:end - 1, :), nodes(2:end, :));
ang = acos(max(min(sum(-v(1:end - 1, :) .* (v(2:end, :)), 2), 1), -1));

newnodes = nodes;
newv = v;
newang = ang;

idx = find(newang < minangle);

while (~isempty(idx))
    newnodes(idx + 1, :) = [];
    newv(idx + 1, :) = [];
    newang(idx) = [];
    idx = unique(idx - (0:(length(idx) - 1))');
    idx1 = idx(idx < size(newnodes, 1));
    newv(idx1, :)  = segvec(newnodes(idx1, :), newnodes(idx1 + 1, :));
    idx1 = idx(idx < size(newv, 1));
    newang(idx1)  = acos(sum(-newv(idx1, :) .* (newv(idx1 + 1, :)), 2));
    idx0 = idx(idx > 1);
    newang(idx0 - 1) = acos(sum(-newv(idx0 - 1, :) .* (newv(idx0, :)), 2));
    idx = find(newang < minangle);
end

if (nargout > 1)
    len = newnodes(1:end - 1, :) - newnodes(2:end, :);
    len = sqrt(sum(len .* len, 2));
end

function v = segvec(n1, n2)

v = n2 - n1;
normals = sqrt(sum(v .* v, 2));
v = v ./ repmat(normals, 1, size(v, 2));
