function [node,face,yz0,yz1]=extrudecurve(xy, yz, Nx, Nz, Nextrap, spacing, anchor, dotopbottom)
%
% [node,face,yz0,yz1]=extrudecurve(xy, yz, Nx, Nz, Nextrap, spacing, anchor)
% 
% create a triangular surface mesh by swining a 2D spline along another 2D
% spline curve
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%      xy: a 2D spline path, along which the surface is extruded, defined
%          on the x-y plane
%      yz: a 2D spline which will move along the path to form a surface,
%          defined on the y-z plane
%      Nx: the count of sample points along the extrusion path (xy), if
%          ignored, it is 40
%      Nz: the count of sample points along the curve to be extruded (yz),
%          if ignored, it is 40
%      Nextrap: number of points to extrapolate outside of the xy/yz
%          curves, 0 if ignored
%      spacing: define a spacing scaling factor for spline interpolations,
%          1 if ignored
%      anchor: the 3D point in the extruded curve plane (yz) that is aligned
%          at the nodes long the extrusion path. this point does not have
%          to be located on the yz curve; orig = [ox oy oz], if ignored, it
%          is set as the point on the interpolated yz with the largested
%          y-value
%      dotopbottom: a flag, if set to 1, tessellated top and bottom faces
%          will be added. default is 0.
%
% output:
%      node: 3D node coordinates for the generated surface mesh
%      face: triangular face patches of the generated surface mesh, each
%           row represents a triangle denoted by the indices of the 3 nodes
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin<3)
    Nx=30;
end
if(nargin<4)
    Nz=30;
end
if(nargin<5)
    Nextrap=0;
end
if(nargin<6)
    spacing=1;
end

xrange=max(xy(:,1))-min(xy(:,1));
dx=xrange/Nx;
xi=min(xy(:,1))-Nextrap*dx:spacing*dx:max(xy(:,1))+Nextrap*dx;
pxy = spline(xy(:,1), xy(:,2));

yi = ppval(pxy,xi);
dy = gradient(yi);
dx = gradient(xi);

nn=sqrt(dx.*dx + dy.*dy);
normaldir=[dx(:)./nn(:) dy(:)./nn(:)];

zrange=max(yz(:,2))-min(yz(:,2));
dz=zrange/Nz;
zi=min(yz(:,2))-Nextrap*dz:spacing*dz:max(yz(:,2))+Nextrap*dz;
pyz = spline(yz(:,2), yz(:,1));

yyi = ppval(pyz,zi);

if(~exist('anchor','var') || isempty(anchor))
    [ymax, loc]=max(yyi);
    anchor=[0 yyi(loc) zi(loc)];
end
    
node=zeros(length(zi)*length(xi),3);
face=zeros(2*(length(zi)-1)*(length(xi)-1),3);

xyz=[yyi(:) yyi(:) zi(:)];
xyz(:,1)=0;
for i=1:length(xi)
    rot2d=[normaldir(i,1) -normaldir(i,2); normaldir(i,2) normaldir(i,1)];
    offset=[xi(i) yi(i) anchor(2)];
    newyz=[((rot2d*(xyz(:,1:2)-repmat(anchor(1:2),size(xyz,1),1))' + repmat(offset(:,1:2),length(zi),1)'))' , xyz(:,end)];
    node((i-1)*length(zi)+1:i*length(zi),:)=newyz;
    if(i>1)
       face((i-2)*2*(length(zi)-1)+1:(i-1)*2*(length(zi)-1),:)= ...
                       [ [1:length(zi)-1 ;  (1-length(zi):-1) ; 2:length(zi)]+(i-1)*length(zi)  ...
                       [ (1-length(zi):-1); (2-length(zi):0)  ; 2:length(zi)]+(i-1)*length(zi) ]';
    end
    if(i==Nextrap+1)
        yz0=newyz(Nextrap+1:end-Nextrap,:);
    end
    if(i==length(xi)-Nextrap)
        yz1=newyz(Nextrap+1:end-Nextrap,:);
    end
end

% add two flat polygons on the top and bottom of the contours 
% to ensure the enclosed surface is not truncated by meshfix

if(nargin>=8 && dotopbottom==1)
    nump = length(xi);
    C = [(1:(nump-1))' (2:nump)'; nump 1];
    dt = delaunayTriangulation(xi(:), yi(:), C);
    io = dt.isInterior();
    endface=dt(io,:);
    endface=(endface-1)*size(newyz,1)+1;

    % append the top/bottom faces to the extruded mesh
    face=[face; endface; endface+size(newyz,1)-1];
end

[node,face]=meshcheckrepair(node,face,'dup');