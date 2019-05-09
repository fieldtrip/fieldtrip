function imgdiff=imedge3d(binimg,isdiff)
%
% imgdiff=imedge3d(binimg,isdiff)
%
% Extract the boundary voxels from a binary image
% 
% author: Aslak Grinsted <ag at glaciology.net>
% modified by Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input: 
%   binimg: a 3D binary image
%   isdiff: if isdiff=1, output will be all voxels which 
%         is different from its neighbors; if isdiff=0 or 
%         ignored, output will be the edge voxels of the 
%         non-zero regions in binimg
%
% output:
%   imgdiff: a 3D logical array with the same size as binimg
%            with 1 for voxels on the boundary and 0 otherwise 
% 
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

invol=1;
if(nargin==2)
	invol=isdiff;
end
binimg=logical(binimg);
imgdiff=xor(binimg,binimg(:,:,[1 1:end-1]));
imgdiff=imgdiff|xor(binimg,binimg(:,:,[2:end end]));
imgdiff=imgdiff|xor(binimg,binimg(:,[1 1:end-1],:));
imgdiff=imgdiff|xor(binimg,binimg(:,[2:end end],:));
imgdiff=imgdiff|xor(binimg,binimg([1 1:end-1],:,:));
imgdiff=imgdiff|xor(binimg,binimg([2:end end],:,:));
if(invol)
	imgdiff=imgdiff&binimg;
end
