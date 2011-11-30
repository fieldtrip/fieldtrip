function [no,el,regions,holes]=v2s(img,isovalues,opt,method)
%
% [no,el,regions,holes]=v2s(img,isovalues,opt,method)
%
% surface mesh generation from binary or gray-scale volumetric images
% shortcut for vol2surf
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% inputs and outputs are similar to those defined in vol2surf
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin==3)
   method='cgalsurf';
end

if(strcmp(method,'cgalmesh'))
   [no,tet,el]=v2m(uint8(img),isovalues,opt,1000,method);
   el=unique(el(:,1:3),'rows');
   [no,el]=removeisolatednode(no(:,1:3),el(:,1:3));
   regions=[]; holes=[];
   return;
end

[no,el,regions,holes]=vol2surf(img,1:size(img,1),1:size(img,2),1:size(img,3),opt,1,method,isovalues);
