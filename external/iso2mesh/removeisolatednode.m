function [no,el]=removeisolatednode(node,elem)
%
% [no,el]=removeisolatednode(node,elem)
%
% remove isolated nodes: nodes that are not included in any element
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%     node: list of node coordinates
%     elem: list of elements of the mesh
%
% output:
%     no: node coordinates after removing the isolated nodes
%     el: element list of the resulting mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

oid=1:size(node,1);       % old node index
idx=setdiff(oid,elem(:)); % indices to the isolated nodes
idx=sort(idx);
delta=zeros(size(oid));   
delta(idx)=1;
delta=-cumsum(delta);     % calculate the new node index after removing the isolated nodes
oid=oid+delta;            % map to new index
el=oid(elem);             % element list in the new index
no=node;                  
no(idx,:)=[];             % remove the isolated nodes
