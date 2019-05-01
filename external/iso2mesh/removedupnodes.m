function [newnode,newelem]=removedupnodes(node,elem,tol)
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

if(nargin>=3 && tol~=0)
    node=round(node/tol)*tol;
end
[newnode,I,J]=unique(node,'rows');
if(iscell(elem))
    newelem=cellfun(@(x) J(x)', elem,'UniformOutput',false);
else
    newelem=J(elem);
end
