function [node,elem]=binsurface(img,nface)
%
% [node,elem]=binsurface(img,nface)
%
% fast isosurface extraction from 3D binary images
%
% author: Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input: 
%   img:  a 3D binary image
%   nface: nface=3 or ignored - for triangular faces, 
%          nface=4 - square faces
%          nface=0 - return a boundary mask image via node
%
% output:
%   elem: integer array with dimensions of NE x nface, each row represents
%         a surface mesh face element 
%   node: node coordinates, 3 columns for x, y and z respectively
%
% the outputs of this subroutine can be easily plotted using 
%     patch('Vertices',node,'faces',elem,'FaceVertexCData',node(:,3),
%           'FaceColor','interp');
% if the surface mesh has triangular faces, one can plot it with
%     trisurf(elem,node(:,1),node(:,2),node(:,3))
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
% 

dim=size(img);
newdim=dim+[1 1 1];

% find the jumps (0->1 or 1->0) for all directions
d1=diff(img,1,1);
d2=diff(img,1,2);
d3=diff(img,1,3);
[ix,iy]=find(d1==1|d1==-1);
[jx,jy]=find(d2==1|d2==-1);
[kx,ky]=find(d3==1|d3==-1);

% compensate the dim. reduction due to diff, and 
% wrap the indices in a bigger array (newdim)
ix=ix+1;
[iy,iz]=ind2sub(dim(2:3),iy);
iy=sub2ind(newdim(2:3),iy,iz);

[jy,jz]=ind2sub([dim(2)-1,dim(3)],jy);
jy=jy+1;
jy=sub2ind(newdim(2:3),jy,jz);

[ky,kz]=ind2sub([dim(2),dim(3)-1],ky);
kz=kz+1;
ky=sub2ind(newdim(2:3),ky,kz);

id1=sub2ind(newdim,ix,iy);
id2=sub2ind(newdim,jx,jy);
id3=sub2ind(newdim,kx,ky);

if(nargin==2 && nface==0)
	elem=[id1 id2 id3];
	node=zeros(newdim);
	node(elem)=1;
	node=node(2:end-1,2:end-1,2:end-1)-1;
	return
end

% populate all the triangles
xy=newdim(1)*newdim(2);

if(nargin==1 || nface==3)  % create triangular elements
    elem=[id1 id1+newdim(1) id1+newdim(1)+xy; id1 id1+newdim(1)+xy id1+xy];
    elem=[elem; id2 id2+1 id2+1+xy; id2 id2+1+xy id2+xy];
    elem=[elem; id3 id3+1 id3+1+newdim(1); id3 id3+1+newdim(1) id3+newdim(1)];
else                       % create box elements
    elem=[id1 id1+newdim(1) id1+newdim(1)+xy id1+xy];
    elem=[elem; id2 id2+1 id2+1+xy id2+xy];
    elem=[elem; id3 id3+1 id3+1+newdim(1) id3+newdim(1)];
end
% compress the node indices
nodemap=zeros(max(elem(:)),1);
nodemap(elem(:))=1;
id=find(nodemap);

nodemap=0;
nodemap(id)=1:length(id);
elem=nodemap(elem);

% create the coordiniates
[xi,yi,zi]=ind2sub(newdim,id);

% assuming the origin [0 0 0] is located at the lower-bottom corner of the image
node=[xi(:),yi(:),zi(:)]-1;
