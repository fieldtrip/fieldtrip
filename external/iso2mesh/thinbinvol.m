function vol=thinbinvol(vol,layer,nobd)
%
% vol=thinbinvol(vol,layer,nobd)
%
% thinning a binary volume by a given pixel width
% this is similar to bwmorph(vol,'thin',n) except 
% this does it in 3d and only run thinning for 
% non-zero elements (and hopefully faster)
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%     vol: a volumetric binary image
%     layer: number of iterations for the thickenining
%     nobd: (optional) if set to 1, boundaries will not 
%            erode. if not given, nobd=0.
%
% output:
%     vol: the volume image after the thinning operations
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

dim=size(vol);
dxy=dim(1)*dim(2);
fulllen=prod(dim);

if(nargin<3)
    nobd=0;
end

if(nobd==1)
    bdmask=vol;
    if(ndims(vol)==2)
        bdmask(2:end-1,2:end-1)=0;
    elseif(ndims(vol)==3)
        bdmask(2:end-1,2:end-1,2:end-1)=0;
    end
end

for i=1:layer
	idx=find(~vol);
	idxnew=[idx+1; idx-1;idx+dim(1);idx-dim(1);idx+dxy;idx-dxy];
	idxnew=idxnew(find(idxnew>0 & idxnew<fulllen));
	vol(idxnew)=0;
        if(nobd)
             vol = vol | bdmask;
        end
end
