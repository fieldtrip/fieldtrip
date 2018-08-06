function [mask, weight]=mesh2mask(node,face,xi,yi,hf)
%
% [mask weight]=mesh2mask(node,face,Nxy)
%   or
% [mask weight]=mesh2mask(node,face,[Nx,Ny])
%   or
% [mask weight]=mesh2mask(node,face,xi,yi,hf)
%
% fast rasterization of a 2D mesh to an image with triangle index labels
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
% date for initial version: July 18,2013
%
% input:
%      node: node coordinates, dimension N by 2 or N by 3 array
%      face: a triangle surface, N by 3 or N by 4 array
%      Nx,Ny,Nxy: output image in x/y dimensions, or both
%      xi,yi: linear vectors for the output pixel center positions in x/y
%      hf: (optional) the handle of a pre-created figure window, for faster 
%          rendering
%
% output:
%      mask: a 2D image, the value of each pixel is the index of the
%            enclosing triangle, if the pixel is outside of the mesh, NaN
%      weight: (optional) a 3 by Nx by Ny array, where Nx/Ny are the dimensions for
%            the mask
%
% note: This function only works in MATLAB when the DISPLAY is not 
%       disabled. The maximum size of the mask output is limited by the 
%       screen size.
%
% example:
%
%   [no,fc]=meshgrid6(0:5,0:5);
%   [mask weight]=mesh2mask(no,fc,-1:0.1:5,0:0.1:5);
%   imagesc(mask);
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin==3 && length(xi)==1 && xi>0)
    mn=min(node);
    mx=max(node);
    df=(mx(1:2)-mn(1:2))/xi;
elseif(nargin==3 && length(xi)==2 && all(xi>0))
    mn=min(node);
    mx=max(node);
    df=(mx(1:2)-mn(1:2))./xi;
elseif(nargin==4 || nargin==5)
    mx=[max(xi) max(yi)];
    mn=[min(xi) min(yi)];
    df=[min(diff(xi(:))) min(diff(yi(:)))];
else
    error('you must give at least xi input');
end
if(size(node,2)<=1 || size(face,2)<=2)
    error('node must have 2 or 3 columns; face can not have less than 2 columns');
end

if(nargin<5)
    hf=figure('visible','on');
else
    clf(hf);
end
patch('Vertices',node,'Faces',face,'linestyle','none','FaceColor','flat',...
 'FaceVertexCData',(1:size(face,1))','CDataMapping', 'scaled');
set(gca, 'Position', [0 0 1 1]);
cm=jet(size(face,1));
colormap(cm);
axis off
set(gca,'xlim',[mn(1) mx(1)]);
set(gca,'ylim',[mn(2) mx(2)]);
set(gca,'clim',[1 size(face,1)]);

output_size = round((mx(1:2)-mn(1:2))./df);%Size in pixels

if(isoctavemesh || isempty(getenv('DISPLAY')))
    resolution = 300; %Resolution in DPI
    set(gcf,'PaperPositionMode','manual')
    set(gcf,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
    deletemeshfile(mwpath('post_mesh2mask.png'));
    print(mwpath('post_mesh2mask.png'),'-dpng',['-r' num2str(resolution)]);
    mask=imread(mwpath('post_mesh2mask.png'));
else
    pos=get(hf,'position');
    pos(3:4)=max(pos(3:4),output_size+20);
    set(hf,'position',pos);
    set(gca, 'Units','pixels','position',[1, 1, output_size(1), output_size(2)]);
    mask=getframe(gca);
    if(any(size(mask.cdata)<[output_size([2 1]) 3]))
        error('the requested rasterization grid is larger than the screen resolution');
    end
    mask=mask.cdata(1:output_size(2),1:output_size(1),:);
end
if(nargin<5)
    close(hf);
end
mask=int32(reshape(mask,[size(mask,1)*size(mask,2) size(mask,3)]));
[isfound,locb]=ismember(mask,floor(cm*255),'rows');
locb(isfound==0)=nan;

mask=rot90(reshape(locb,output_size([2 1]))');

if(nargout>=2)
    xi=mn(1)+df(1)/2:df(1):mx(1);
    yi=mn(2)+df(2)/2:df(2):mx(2);
    weight=barycentricgrid(node,face,xi,yi,mask);
    if(size(face,2)>=4)
        badidx=find(weight(1,:,:)<0 | weight(2,:,:)<0 | weight(3,:,:)<0);
        badidx=badidx(face(mask(badidx),3)~=face(mask(badidx),4));
        weight2=barycentricgrid(node,face(:,[1 3 4]),xi,yi,mask);
        weight(:,badidx)=0;
        weight([1 3 4],badidx)=weight2(:,badidx);
    end
end

function weight=barycentricgrid(node,face,xi,yi,mask)
[xx,yy]=meshgrid(xi,yi);
idx=find(~isnan(mask));
eid=mask(idx);
t1=node(face(:,1),:);
t2=node(face(:,2),:);
t3=node(face(:,3),:);
tt=(t2(:,2)-t3(:,2)).*(t1(:,1)-t3(:,1))+(t3(:,1)-t2(:,1)).*(t1(:,2)-t3(:,2));
w(:,1)=(t2(eid,2)-t3(eid,2)).*(xx(idx)-t3(eid,1))+(t3(eid,1)-t2(eid,1)).*(yy(idx)-t3(eid,2));
w(:,2)=(t3(eid,2)-t1(eid,2)).*(xx(idx)-t3(eid,1))+(t1(eid,1)-t3(eid,1)).*(yy(idx)-t3(eid,2));
w(:,1)=w(:,1)./tt(eid);
w(:,2)=w(:,2)./tt(eid);
w(:,3)=1-w(:,1)-w(:,2);
weight=zeros(3,size(mask,1),size(mask,2));
ww=zeros(size(mask));
ww(idx)=w(:,1);
weight(1,:,:)=ww;
ww(idx)=w(:,2);
weight(2,:,:)=ww;
ww(idx)=w(:,3);
weight(3,:,:)=ww;
