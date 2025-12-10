function [p, e0] = ray2surf(node, elem, p0, v0, e0)
%
% [p,e0]=ray2surf(node,elem,p0,v0,e0)
%
% Determine the entry position and element for a ray to intersect a mesh
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     node: the mesh coordinate list
%     elem: the tetrahedral mesh element list, 4 columns
%     p0: origin of the ray
%     v0: direction vector of the ray
%     e0: search direction: '>' forward search, '<' backward search, '-' bidirectional
%
% output:
%     p: the intersection position
%     e0: if found, the index of the intersecting element ID
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

p = p0;
if (size(elem, 2) == 3)
    face = elem;
else
    face = volface(elem);
end
[t, u, v, idx] = raytrace(p0, v0, node, face);
if (isempty(idx))
    error('ray does not intersect with the mesh');
else
    t = t(idx);
    if (e0 == '>')
        %         idx1=find(t>=0);
        idx1 = find(t >= 1e-10);
    elseif (e0 == '<')
        idx1 = find(t <= 0);
    elseif (isnan(e0) || e0 == '-')
        idx1 = 1:length(t);
    else
        error('ray direction specifier is not recognized');
    end
    if (isempty(idx1))
        error('no intersection is found along the ray direction');
    end
    t0 = abs(t(idx1));
    [tmin, loc] = min(t0);
    faceidx = idx(idx1(loc));

    % update source position
    p = p0 + t(idx1(loc)) * v0;

    if (nargout < 2)
        return
    end

    % find initial element id
    if (size(elem, 2) == 3)
        e0 = faceidx;
    else
        felem = sort(face(faceidx, :));
        f = elem;
        f = [f(:, [1, 2, 3])
             f(:, [2, 1, 4])
             f(:, [1, 3, 4])
             f(:, [2, 4, 3])];
        [tf, loc] = ismember(felem, sort(f, 2), 'rows');
        loc = mod(loc, size(elem, 1));
        if (loc == 0)
            loc = size(elem, 1);
        end
        e0 = loc;
    end
end
