function [faces, idx, facemap] = uniqfaces(elem)
%
% [faces,idx,facemap]=uniqfaces(elem)
%
% return the unique face list from a or tetrahedral mesh
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%     elem: a list of elements, each row is a list of nodes for an element.
%           elem can have 2, 3 or 4 columns
%
% output:
%     face: unique faces in the mesh, denoted by a triplet of node indices
%     idx:  index of the output in the raw face list (returned by meshface)
%     facemap: index of the raw faces in the output list (for triangular mesh)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if (size(elem) == 3)
    faces = elem;
elseif (size(elem) >= 4)
    faces = meshface(elem);
else
    error('invalid input');
end

[ufaces, idx, jdx] = unique(sort(faces, 2), 'rows');
faces = faces(idx, :);
if (nargout > 2)
    facemap = reshape(jdx, [size(elem, 1) nchoosek(size(elem, 2), 3)]);
end
