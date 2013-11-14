function vol=thinbinvol(vol,layer)
%
% vol=thinbinvol(vol,layer)
%
% thinning a binary volume by a given pixel width
% this is similar to bwmorph(vol,'thin',n) except 
% this does it in 3d and only does thinning for 
% non-zero elements (and hopefully faster)
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%     vol: a volumetric binary image
%     layer: number of iterations for the thickenining
%
% output:
%     vol: the volume image after the thickening
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

dim=size(vol);
dxy=dim(1)*dim(2);
fulllen=prod(dim);

for i=1:layer
	idx=find(~vol);
	idxnew=[idx+1; idx-1;idx+dim(1);idx-dim(1);idx+dxy;idx-dxy];
	idxnew=idxnew(find(idxnew>0 & idxnew<fulllen));
	vol(idxnew)=0;
end
