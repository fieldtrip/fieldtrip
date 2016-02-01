function [newnode,newelem]=removedupnodes(node,elem)
%
% [newnode,newelem]=removedupnodes(node,elem)
%
% removing the duplicated nodes from a mesh
%
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input:
%   elem: integer array with dimensions of NE x 4, each row contains
%         the indices of all the nodes for each tetrahedron
%   node: node coordinates, 3 columns for x, y and z respectively
%
% output:
%   newnode: nodes without duplicates
%   newelem: elements with only the unique nodes
%
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[newnode,I,J]=unique(node,'rows');
newelem=J(elem);
