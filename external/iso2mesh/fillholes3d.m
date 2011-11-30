function resimg=fillholes3d(img,ballsize)
%
% resimg=fillholes3d(img,ballsize)
%
% close a 3D image with the speicified gap size and then fill the holes
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% 
% input:
%    img: a 3D binary image
%    ballsize: maximum gap size for image closing
%
% output:
%    resimg: the image free of holes
%
% this function requires the image processing toolbox for matlab/octave
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(ballsize)
  resimg = imclose(img,strel(ones(ballsize,ballsize,ballsize)));
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
