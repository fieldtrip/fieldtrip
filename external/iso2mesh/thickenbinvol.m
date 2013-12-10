function vol=thickenbinvol(vol,layer)
%
% vol=thickenbinvol(vol,layer)
%
% thickening a binary volume by a given pixel width
% this is similar to bwmorph(vol,'thicken',3) except 
% this does it in 3d and only does thickening for 
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

offs=[1,-1,dim(1),-dim(1),dxy,-dxy];
for i=1:layer
	idx=find(vol);
    for j=1:6
    	idxnew=idx+offs(j);
        idxnew=idxnew(find(idxnew>0 & idxnew<fulllen));
    	vol(idxnew)=1;
    end
end
