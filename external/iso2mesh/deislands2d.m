function cleanimg=deislands2d(img,sizelim)
%
% cleanimg=deislands2d(img,sizelim)
%
% remove isolated islands on a 2D image below speicified size limit
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
% 	img: a 2D binary image
% 	sizelim: a integer as the maximum pixel size of a isolated region
%
% output:
% 	cleanimg: a binary image after removing islands below sizelim
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

img=squeeze(img);
maxisland=-1;
if(nargin==2) maxisland=sizelim; end

islands={};

cleanimg=zeros(size(img));
if(sum(img(:)))
    img=imclose(img, strel('disk',3));
    islands=bwislands(img);
end

if(length(islands))
    % remove small islands of the foreground
    maxblock=-1;
    maxblockid=-1;
    if(maxisland<0)
      for i=1:length(islands)
        if(length(islands{i})>maxblock)
            maxblockid=i;
            maxblock=length(islands{i});
        end
      end
      if(maxblock>0)
          cleanimg(islands{maxblockid})=1;
      end
    else
      for i=1:length(islands)
        if(length(islands{i})>maxisland)
            cleanimg(islands{i})=1;
        end
      end
    end

    % remote small islands of the background
    
    cleanimg=imfill(cleanimg,'holes');
end
