function [newpt elemid weight] = proj2mesh(v, f, pt, nv, cn, radmax)
%  [newpt elemid weight]=proj2mesh(v,f,pt,nv,cn)
%
%  project a point cloud on to the surface mesh (surface can only be triangular)
%
%  author: Qianqian Fang <q.fang at neu.edu>
%  date: 12/12/2008
%
% parameters:
%      v: node coordinate of the surface mesh (nn x 3)
%      f: element list of the surface mesh (3 columns for
%            triangular mesh, 4 columns for cubic surface mesh)
%      pt: points to be projected, 3 columns for x,y and z respectively
%      nv: nodal norms (vector) calculated from nodesurfnorm.m
%          with dimensions of (size(v,1),3)
%      cn: a integer vector with the length of p, denoting the closest
%          surface nodes (indices of v) for each point in p. this
%          value can be calculated from dist2surf.m
%      radmax: if speicified, the search for elements to project will be
%          limited to those within a bounding box with half-edge-length
%          of radmax centered at the point to be projected
%
%      if nv and cn are not supplied, proj2mesh will project the point
%      cloud onto the surface by the direction pointing to the centroid
%      of the mesh
%
% outputs:
%      newpt: the projected points from p
%      elemid: a vector of length of p, denotes which surface trangle (in elem)
%             contains the projected point
%      weight: the barycentric coordinate for each projected points, these are
%             the weights
%
% Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%
% this function is part of "metch" toobox, see COPYING for license

cent = mean(v);
enum = length(f);
ec = reshape(v(f(:, 1:3)', :)', [3 3, enum]);
centroid = squeeze(mean(ec, 2));
newpt = zeros(size(pt, 1), 3);
elemid = zeros(size(pt, 1), 1);
weight = zeros(size(pt, 1), 3);
[idoldmesh, loc] = ismember(pt, v, 'rows');
idnode = find(idoldmesh);
if (~isempty(idnode))
    [tt, ll] = ismember(loc(idnode), f);
    [p1, p2] = ind2sub(size(f), ll); % p1 is the index in f
    newpt(idnode, :) = pt(idnode, :);
    elemid(idnode) = p1;
    weight(sub2ind(size(weight), idnode, p2)) = 1;
end
radlimit = -1;

if (nargin >= 5)
    % if nv and cn are supplied, use nodal norms to project the points
    direction = nv(cn, :);
    if (nargin >= 6)
        radlimit = radmax;
    end
elseif (nargin == 3)
    % otherwise, project toward the centroid
    direction = pt - repmat(cent, size(pt, 1), 1);
end

for t = 1:size(pt, 1)
    if (idoldmesh(t) ~= 0)
        continue
    end
    maxdist = sqrt(sum((pt(t, :) - cent) .* (pt(t, :) - cent)));

    if (radlimit > 0)
        maxdist = radlimit;
    end

    idx = find(centroid(1, :) > pt(t, 1) - maxdist & ...
               centroid(1, :) < pt(t, 1) + maxdist & ...
               centroid(2, :) > pt(t, 2) - maxdist & ...
               centroid(2, :) < pt(t, 2) + maxdist & ...
               centroid(3, :) > pt(t, 3) - maxdist & ...
               centroid(3, :) < pt(t, 3) + maxdist);

    dist = centroid(:, idx);
    dist(1, :) = dist(1, :) - pt(t, 1);
    dist(2, :) = dist(2, :) - pt(t, 2);
    dist(3, :) = dist(3, :) - pt(t, 3);
    c0 = sum(dist .* dist);

    % sort the distances to accelate the calculation
    [c1, sorted] = sort(c0);

    for i = 1:length(idx)

        % project the point along vector direction and calculate the intersection to a plane

        [inside, p, w] = linextriangle(pt(t, :), pt(t, :) + direction(t, :), v(f(idx(sorted(i)), :), :));

        % the intersection is within the current trianglar surface
        if (inside)
            newpt(t, :) = p;
            weight(t, :) = w;
            elemid(t) = idx(sorted(i));
            break
        end
    end
end
