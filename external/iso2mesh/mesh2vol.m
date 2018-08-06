function [mask weight]=mesh2vol(node,elem,xi,yi,zi)
%
% [mask weight]=mesh2vol(node,face,Nxyz)
% [mask weight]=mesh2vol(node,face,[Nx,Ny,Nz])
% [mask weight]=mesh2vol(node,face,xi,yi,zi,hf)
%   or
% newval=mesh2vol(node_val,face,...)
%
% fast rasterization of a 3D mesh to a volume with tetrahedron index labels
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
% date for initial version: Feb 10,2014
%
% input:
%      node: node coordinates, dimension N by 2 or N by 3 array
%      nodeval: a 4-column array, the first 3 columns are the node coordinates, 
%            the last column denotes the values associated with each node
%      face: a triangle surface, N by 3 or N by 4 array
%      Nx,Ny,Nxy: output image in x/y dimensions, or both
%      xi,yi: linear vectors for the output pixel center positions in x/y
%      hf: the handle of a pre-created figure window for faster rendering
%
% output:
%      mask: a 3D image, the value of each pixel is the index of the
%            enclosing triangle, if the pixel is outside of the mesh, NaN
%      weight: (optional) a 3 by Nx by Ny array, where Nx/Ny are the dimensions
%            for the mask
%      newval: when the node has 4 columns, the last column represents the
%            values (color) at each node, the output newval is the rasterized
%            mesh value map over the specified grid.
%
% note: This function only works for matlab
%
% example:
%
%   [no,el]=meshgrid6(0:5,0:5,0:3);
%   mask=mesh2vol(no,el(:,1:4),0:0.1:5,0:0.1:5,0:0.1:4);
%   imagesc(mask(:,:,3))
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

nodeval=[];
if(size(node,2)==4)
   nodeval=node(:,4);
   node(:,4)=[];
end

if(nargin==3 && length(xi)==1 && xi>0)
    mn=min(node);
    mx=max(node);
    df=(mx(1:min(3,size(node,2)))-mn(1:min(3,size(node,2))))/xi;
elseif(nargin==3 && length(xi)==3 && all(xi>0))
    mn=min(node);
    mx=max(node);
    df=(mx(1:min(3,size(node,2)))-mn(1:min(3,size(node,2))))./(xi(:)');
elseif(nargin==5)
    mx=[max(xi) max(yi) max(zi)];
    mn=[min(xi) min(yi) min(zi)];
    df=[min(diff(xi(:))) min(diff(yi(:))) min(diff(zi(:)))];
else
    error('you must give at least xi input');
end

xi=mn(1):df(1):mx(1);
yi=mn(2):df(2):mx(2);
zi=mn(3):df(3):mx(3);

if(size(node,2)~=3 || size(elem,2)<=3)
    error('node must have 3 columns; face can not have less than 3 columns');
end

nz=length(zi);
mask=zeros(length(xi)-1,length(yi)-1,length(zi)-1);
if(nargout>1 || ~isempty(nodeval))
   weight=zeros(4,length(xi)-1,length(yi)-1,length(zi)-1);
end

hf=figure('visible','on');
for i=1:nz
    if(~isempty(nodeval))
        [cutpos,cutvalue,facedata,elemid]=qmeshcut(elem,node,nodeval,['z=' num2str(zi(i))]);
    else
        [cutpos,cutvalue,facedata,elemid]=qmeshcut(elem,node,node(:,1),['z=' num2str(zi(i))]);
    end
    if(isempty(cutpos))
        continue;
    end
    if(nargout>1 || ~isempty(nodeval))
        [maskz, weightz]=mesh2mask(cutpos,facedata,xi,yi,hf);
        weight(:,:,:,i)=weightz;
    else
        maskz=mesh2mask(cutpos,facedata,xi,yi,hf);
    end
    idx=find(~isnan(maskz));
    if(~isempty(nodeval))
        eid=facedata(maskz(idx),:);
        maskz(idx)=cutvalue(eid(:,1)).*weightz(1,idx)'+cutvalue(eid(:,2)).*weightz(2,idx)'+...
            cutvalue(eid(:,3)).*weightz(3,idx)'+cutvalue(eid(:,4)).*weightz(4,idx)';
    else
        maskz(idx)=elemid(maskz(idx));
    end
    mask(:,:,i)=maskz;
end
close(hf);
