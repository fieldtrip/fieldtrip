function quality=meshquality(node,elem,maxnode)
%
% quality=meshquality(node,elem)
%
% compute the Joe-Liu mesh quality measure of an N-D mesh (N<=3)
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2011/02/26
%
% input:
%    node:  node coordinates of the mesh (nn x 3)
%    elem:  element table of an N-D mesh (ne x (N+1))
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
%                    BIT 34 (2) (1994) 268-287.
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin<3)
    maxnode=4;
end
if(size(elem,2)>maxnode)
    elem=elem(:,1:maxnode);
end
enum=size(elem,1);
vol=elemvolume(node,elem);
edges=meshedge(elem);
ed=node(edges(:,1),:)-node(edges(:,2),:);
ed=sum((ed.*ed)');
ed=sum(reshape(ed,[enum length(ed)/enum])')';
dim=size(elem,2)-1;

coeff=10/9; % for tetrahedral
if(dim==2)
    coeff=1;
end
quality=coeff*dim*2^(2*(1-1./dim))*3^((dim-1)/2)*vol.^(2/dim)./ed;
maxquality=max(quality);
if(maxquality>1)
    quality=quality./maxquality;
end
