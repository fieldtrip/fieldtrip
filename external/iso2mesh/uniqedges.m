function [edges,idx,edgemap]=uniqedges(elem)
%
% [edges,idx,edgemap]=uniqedges(elem)
%
% return the unique edge list from a surface or tetrahedral mesh
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%     elem: a list of elements, each row is a list of nodes for an element.
%           elem can have 2, 3 or 4 columns
%
% output:
%     edge: unique edges in the mesh, denoted by a pair of node indices
%     idx:  index of the output in the raw edge list (returned by meshedge)
%     edgemap: index of the raw edges in the output list (for triangular mesh)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(size(elem)==2)
   edges=elem;
elseif(size(elem)>=3)
   edges=meshedge(elem);
else
   error('invalid input');
end

[uedges,idx,jdx]=unique(sort(edges,2),'rows');
edges=edges(idx,:);
edgemap=reshape(jdx,size(elem));
