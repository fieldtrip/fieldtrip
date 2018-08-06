function [node,face,elem]=meshacylinder(c0,c1,r,tsize,maxvol,ndiv)
%
% [node,face]=meshacylinder(c0,c1,r,tsize,maxvol,ndiv)
%    or
% [node,face,elem]=meshacylinder(c0,c1,r,tsize,maxvol,ndiv)
%
% create the surface and (optionally) tetrahedral mesh of a 3D cylinder
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input: 
%   c0, c1:  cylinder axis end points
%   r:   radius of the cylinder
%   tsize: maximum surface triangle size on the sphere
%   maxvol: maximu volume of the tetrahedral elements
%   ndiv: approximate the cylinder surface into ndiv flat pieces, if 
%         ignored, ndiv is set to 20
%
% output:
%   node: node coordinates, 3 columns for x, y and z respectively
%   face: integer array with dimensions of NB x 3, each row represents
%         a surface mesh triangle 
%   elem: (optional) integer array with dimensions of NE x 4, each row 
%         represents a tetrahedron 
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
% 

if(nargin<3)
    error('you must at least provide c0, c1, and r');
end
if(r<=0 || all(c0==c1))
    error('invalid cylinder parameters');
end
c0=c0(:);
c1=c1(:);
len=sqrt(sum((c0-c1).*(c0-c1)));
% define the axial vector v0 and a perpendicular vector t
v0=c1-c0;

% calculate the cylinder end face nodes
if(nargin<6) ndiv=20; end

dt=2*pi/ndiv;
theta=dt:dt:2*pi;
cx=r*cos(theta);
cy=r*sin(theta);
p0=[cx(:) cy(:) zeros(ndiv,1)];
p1=[cx(:) cy(:) len*ones(ndiv,1)];
pp=[p0;p1];
no=rotatevec3d(pp,v0)+repmat(c0',size(pp,1),1);

count=1;
for i=1:ndiv-1
   fc{count}={[i i+ndiv i+ndiv+1 i+1],1}; count=count+1;
end
i=ndiv;
fc{count}={[i i+ndiv 1+ndiv 1],1}; count=count+1;
fc{count}={1:ndiv,2};count=count+1;  % bottom inner circle
fc{count}={1+ndiv:2*ndiv,3};count=count+1;  % top inner circle

if(nargin==3)
    tsize=len/10;
end
if(nargin<5)
    maxvol=tsize*tsize*tsize;
end
[node,elem,face]=surf2mesh(no,fc,min(no),max(no),1,maxvol,[0 0 1],[],0);
