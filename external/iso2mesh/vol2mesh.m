function [node,elem,face,regions]=vol2mesh(img,ix,iy,iz,opt,maxvol,dofix,method,isovalues)
%
% [node,elem,face,regions]=vol2mesh(img,ix,iy,iz,opt,maxvol,dofix,method,isovalues)
%
% convert a binary (or multi-valued) volume to tetrahedral mesh
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%	 img: a volumetric binary image
%	 ix,iy,iz: subvolume selection indices in x,y,z directions
%	 opt: as defined in vol2surf.m
%	 maxvol: target maximum tetrahedral elem volume
%	 dofix: 1: perform mesh validation&repair, 0: skip repairing
%	 method: 'cgalsurf' or omit: use CGAL surface mesher
%		 'simplify': use binsurface and then simplify
%		 'cgalmesh': use CGAL 3.5 3D mesher for direct mesh generation [new]
%
%		 generally speaking, 'cgalmesh' is the most robust path
%		 if you want to product meshes from binary or multi-region
%		 volumes, however, its limitations include 1) only accept 
%		 uint8 volume, and 2) can not extract meshes from gray-scale
%		 volumes. If ones goal is to process a gray-scale volume,
%		 he/she should use the 'cgalsurf' option. 'simplify' approach
%		 is not recommended unless other options failed.
%	 isovalues: a list of isovalues where the levelset is defined
%
% output:
%	 node: output, node coordinates of the tetrahedral mesh
%	 elem: output, element list of the tetrahedral mesh, the last 
%	       column is the region ID
%	 face: output, mesh surface element list of the tetrahedral mesh
%	       the last column denotes the boundary ID
%    region: optional output. if opt.autoregion is set to 1, region
%          saves the interior points for each closed surface component
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin>=8)
	if(strcmp(method,'cgalmesh'))
		vol=img(ix,iy,iz);
		if(length(unique(vol(:)))>64 & dofix==1)
			error([ 'it appears that you are processing a ' ...
                                'grayscale image. Currently cgalmesher ' ...
                                'does not support grayscale images. ' ...
                                'Please use "cgalsurf" method to mesh a grayscale ' ...
                                'volume. If you are certain to run cgalmesher ' ...
                                'on your data, please set dofix=0 and run this again.' ]);
		end
		[node elem,face]=cgalv2m(vol,opt,maxvol);
		return;
	end
end

%first, convert the binary volume into isosurfaces
if(nargin==8)
	[no,el,regions,holes]=vol2surf(img,ix,iy,iz,opt,dofix,method);
elseif(nargin==9)
	[no,el,regions,holes]=vol2surf(img,ix,iy,iz,opt,dofix,method,isovalues);
else
        [no,el,regions,holes]=vol2surf(img,ix,iy,iz,opt,dofix,'cgalsurf');
end
%then, create volumetric mesh from the surface mesh
if(nargin>=8)
   if(strcmp(method,'cgalpoly'))
	[node,elem,face]=cgals2m(no(:,1:3),el(:,1:3),opt,maxvol);
        return;
   end
end

[node,elem,face]=surf2mesh(no,el,[],[],1,maxvol,regions,holes);
