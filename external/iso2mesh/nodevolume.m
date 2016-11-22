function nodevol=nodevolume(node,elem)
%
% nodevol=nodevolume(node,elem)
%
% calculate the Voronoi volume of each node in a simplex mesh
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2009/12/31
%
% input:
%    node:  node coordinates
%    elem:  element table of a mesh
%
% output:
%    nodevol:   volume values for all nodes
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

dim=4;
if(size(elem,2)==3) dim=3; end

vol=elemvolume(node,elem(:,1:dim));

elemnum=size(elem,1);
nodenum=size(node,1);
nodevol=zeros(nodenum,1);
for i=1:elemnum
      nodevol(elem(i,1:dim))=nodevol(elem(i,1:dim))+vol(i);
end
nodevol=nodevol/dim;
