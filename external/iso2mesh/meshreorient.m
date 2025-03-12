function [newelem, evol, idx] = meshreorient(node, elem)
%
% [newelem, evol, idx]=meshreorient(node,elem)
%
% reorder nodes in a surface or tetrahedral mesh to ensure all
% elements are oriented consistently
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2010/05/05
%
% input:
%    node: list of nodes
%    elem: list of elements (each row are indices of nodes of each element)
%    idx: indices of the elements there had negative volume
%
% output:
%    newelem: the element list with consistent ordering
%    evol: the signed element volume before reorientation
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

% calculate the canonical volume of the element (can be a 2D or 3D)
evol = elemvolume(node, elem, 'signed');

% make sure all elements are positive in volume
idx = find(evol < 0);
elem(idx, [end - 1, end]) = elem(idx, [end, end - 1]);
newelem = elem;
