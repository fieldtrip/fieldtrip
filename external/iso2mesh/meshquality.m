function quality=meshquality(node,elem)
%
% quality=meshquality(node,elem)
%
% compute Joe-Liu mesh quality measure of a tetrahedral mesh
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2011/02/26
%
% input:
%    node:  node coordinates of the mesh (nn x 3)
%    elem:  element table of a tetrahedral mesh (ne x 4)
%
% output
%    edge:  edge list; each row is an edge, specified by the starting and
%           ending node indices, the total edge number is
%           size(elem,1) x nchoosek(size(elem,2),2). All edges are ordered
%           by looping through each element first. 
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(size(elem,2)>4)
    elem=elem(:,1:4);
end
enum=size(elem,1);
vol=elemvolume(node,elem);
edges=meshedge(elem);
ed=node(edges(:,1),:)-node(edges(:,2),:);
ed=sum((ed.*ed)');
ed=sum(reshape(ed,[enum length(ed)/enum])')';

quality=12*((3*vol).^(2/3))./ed;
