function cleanimg=deislands3d(img,sizelim)
%
% cleanimg=deislands3d(img,sizelim)
%
% remove isolated islands for 3D image (for each slice)
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%      img: a 3D volumetric image
%      sizelim: maximum island size (in pixels) for each x/y/z slice
%
% output:
%      cleanimg: 3D image after removing the islands
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

maxisland=-1;
if(nargin==2) maxisland=sizelim; end

for i=1:size(img,1)
    if(mod(i,10)==0) fprintf(1,'processing slice x=%d\n',i); end
    img(i,:,:)=deislands2d(img(i,:,:),maxisland);
end
for i=1:size(img,2)
    if(mod(i,10)==0) fprintf(1,'processing slice y=%d\n',i); end
    img(:,i,:)=deislands2d(img(:,i,:),maxisland);
end
for i=1:size(img,3)
    if(mod(i,10)==0) fprintf(1,'processing slice z=%d\n',i); end
    img(:,:,i)=deislands2d(img(:,:,i),maxisland);
end

cleanimg=img;
