function quality=meshquality(node,elem)
%
% quality=meshquality(node,elem)
%
% compute the Joe-Liu mesh quality measure of a tetrahedral mesh
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2011/02/26
%
% input:
%    node:  node coordinates of the mesh (nn x 3)
%    elem:  element table of a tetrahedral mesh (ne x 4)
%
% output:
%    quality: a vector of the same length as size(elem,1), with 
%           each element being the Joe-Liu mesh quality metric (0-1) of 
%           the corresponding element. A value close to 1 represents
%           higher mesh quality (1 means equilateral tetrahedron); 
%           a value close to 0 means nearly degenerated element.
%
% reference:
%    A. Liu, B. Joe, Relationship between tetrahedron shape measures, 
%                    BIT 34 (2) (1994) 268–287.
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
