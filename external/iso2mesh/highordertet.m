function [newnode,newelem]=highordertet(node,elem,order)
%
% [newnode,newelem]=highordertet(node,elem)
%
% generate high-order straight-edge tetrahedral mesh from
% the 1st order tetrahedral mesh
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    node: list of nodes
%    elem: list of elements (each row are indices of nodes of each element)
%    order: optional, the order of the generated mesh; if missing, order=2
%
% output:
%    newnode: all new edge-nodes on the output mesh
%    newelem: the indices of the edge nodes for each original tet element
%
%    currently, this function only supports order=2
%    to combine the newnode/newelem with the old mesh, one should use
%
%    elemfull=[elem(:,1:4) newelem+size(node,1)]; % 10-node element
%    nodefull=[node;newnode];
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin<3)
    order=2;
end
if(order>=3 || order<=1)
    error('currently this function only supports order=2');
end
[edges,idx,newelem]=uniqedges(elem(:,1:4));
newnode=node(edges',1:3);
newnode=reshape(newnode',[3,2,size(edges,1)]);
newnode=squeeze(mean(permute(newnode,[3 2 1]),2));
