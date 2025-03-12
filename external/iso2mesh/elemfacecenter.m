function [newnode, newelem] = elemfacecenter(node, elem)
%
% [newnode,newelem]=elemfacecenter(node,elem)
%
% generate barycentric dual-mesh face center nodes and indices per element
% very similar to highordertet which finds edge-centers instead of
% face-centers
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    node: list of nodes
%    elem: list of elements (each row are indices of nodes of each element)
%
% output:
%    newnode: all new face-nodes on the output mesh
%    newelem: the indices of the face nodes for each original tet element
%
%    to combine the newnode/newelem with the old mesh, one should use
%
%    elemfull=[elem(:,1:4) newelem+size(node,1)];
%    nodefull=[node;newnode];
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[faces, idx, newelem] = uniqfaces(elem(:, 1:4));
newnode = node(faces', 1:3);
newnode = reshape(newnode', [3, 3, size(faces, 1)]);
newnode = squeeze(mean(permute(newnode, [3 2 1]), 2));
