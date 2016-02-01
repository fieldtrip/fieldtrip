function [node,face,elem]=meshunitsphere(tsize,maxvol)
%
% [node,face,elem]=meshunitsphere(tsize,maxvol)
%
% create the surface and/or volumetric mesh of a unit sphere 
% centered at [0 0 0] and radius 1
%
% author: Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input: 
%   tsize: maximum size of the surface triangles (from 0 to 1)
%   maxvol: maximum volume of the tetrahedron; if one wants to return
%           elem without specifying maxvol, maxvol=tsize^3
%
% output:
%   node: node coordinates, 3 columns for x, y and z respectively
%   face: integer array with dimensions of NB x 3, each row represents
%         a surface mesh face element 
%   elem: integer array with dimensions of NE x 4, each row represents
%         a tetrahedron. If ignored, this function only produces the surface
%
% example:
%   [node,face]=meshunitsphere(0.05);
%   [node,face,elem]=meshunitsphere(0.05,0.01);
%   plotmesh(node,elem,'x>0'); axis equal;
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
% 

dim=60;
esize=tsize*dim;
thresh=dim/2-1;
[xi,yi,zi]=meshgrid(0:0.5:dim,0:0.5:dim,0:0.5:dim);
dist=thresh-sqrt((xi-30).^2+(yi-30).^2+(zi-30).^2);
dist(find(dist<0))=0;
clear xi yi zi;

% extract a level-set at v=thresh, being a sphere with R=thresh
% the maximum element size of the surface triangles is tsize*dim

[node,face]=vol2restrictedtri(dist,1,[dim dim dim],dim*dim*dim,30,esize,esize,40000);
node=(node-0.5)*0.5;

[node,face]=removeisolatednode(node,face);

node=(node-30)/28;
r0=sqrt(sum((node.*node)'));
node=node.*repmat(1./r0(:),1,3);

if(nargout==3)
   if(nargin==1) maxvol=tsize*tsize*tsize; end
   [node,elem,face]=surf2mesh(node,face,[-1 -1 -1]*1.1,[1 1 1]*1.1,1,maxvol,[],[]);
   elem=elem(:,1:4);
end
face=face(:,1:3);
