function [no,el,fc,nodemap]=sortmesh(origin,node,elem,ecol,face,fcol)
%
% [no,el,fc]=sortmesh(origin,node,elem,face)
%
% sort nodes and elements in a mesh so that the indexed
% nodes and elements are closer to each order
% (this may reduce cache-miss in a calculation)
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2010/05/06
%
% input:
%    origin: sorting all nodes and elements with the distance and
%            angles wrt this location, if origin=[], it will be 
%            node(1,:)
%    node: list of nodes
%    elem: list of elements (each row are indices of nodes of each element)
%    ecol: list of columns in elem to participate sorting
%    face: list of surface triangles (this can be omitted)
%    fcol: list of columns in face to participate sorting
%
% output:
%    no: node coordinates in the sorted order
%    el: the element list in the sorted order
%    fc: the surface triangle list in the sorted order (can be ignored)
%    nodemap: the new node mapping order, no=node(nodemap,:)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(isempty(origin))
   origin=node(1,:);
end
sdist=node-repmat(origin,size(node,1),1);
[theta,phi,R]=cart2sph(sdist(:,1),sdist(:,2),sdist(:,3));
sdist=[R,phi,theta];
[nval,nodemap]=sortrows(sdist);
no=node(nodemap,:);

[nval,nidx]=sortrows(nodemap);
el=elem;
if(nargin<4)
   ecol=1:size(elem,2);
end
el(:,ecol)=sort(nidx(elem(:,ecol)),2);
el=sortrows(el,ecol);

if(nargin>=5 && nargout==3)
  if(nargin<6)
     fcol=1:size(face,2);
  end
  fc=face;
  fc(:,fcol)=sort(nidx(face(:,fcol)),2);
  fc=sortrows(fc,fcol);
end
