function resimg=fillholes3d(img,maxgap)
%
% resimg=fillholes3d(img,maxgap)
%
% close a 3D image with the speicified gap size and then fill the holes
%
% author: Qianqian Fang, <q.fang at neu.edu>
% 
% input:
%    img: a 3D binary image
%    maxgap: maximum gap size for image closing
%
% output:
%    resimg: the image free of holes
%
% this function requires the image processing toolbox for matlab/octave
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(maxgap)
  resimg = imclose(img,strel(ones(maxgap,maxgap,maxgap)));
else
  resimg=img;
end

if(isoctavemesh)
  if(~exist('bwfill'))
    error('you need to install octave-image toolbox first');
  end
  for i=1:size(resimg,3)
    resimg(:,:,i)=bwfill(resimg(:,:,i),'holes');
  end
else
  resimg=imfill(resimg,'holes');
end
