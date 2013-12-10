function [node,elem,face]=volmap2mesh(img,ix,iy,iz,elemnum,maxvol,thickness,Amat,Bvec)
%
% [node,elem,face]=volmap2mesh(img,ix,iy,iz,thickness,elemnum,maxvol,A,B)
%
% convert a binary volume to tetrahedral mesh followed by an Affine transform
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
% date:   2008/01/12
%
% input: 
%        img, ix,iy,iz, elemnum and  maxvol: see vol2mesh.m
%        thickness: scale z-dimension of the mesh to specified thickness, 
%                   if thickness==0, scaling is bypassed
%        Amat: a 3x3 transformation matrix
%        Bvec: a 3x1 vector
%        Amat and Bvec maps the image index space to real world coordnate system by
%                   [x,y,z]_new=Amat*[x,y,z]_old+Bvec
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[node,elem,face]=vol2mesh(img,ix,iy,iz,elemnum,maxvol);

node(:,1:3)=(Amat*node(:,1:3)'+repmat(Bvec(:),1,size(node,1)))';

if(thickness)
	zmin=min(node(:,3));
	zmax=max(node(:,3));
	node(:,3)=(node(:,3)-zmin)/(zmax-zmin)*thickness;
end
