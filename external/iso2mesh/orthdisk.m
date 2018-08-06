function node=orthdisk(c0,c1,r,ndiv)
%
% node=orthdisk(c0,c1,r,ndiv)
%
% Defining a 3D disk that is orthogonal to the vector c1-c0 
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%     c0: a 1x3 vector for the origin
%     c1: a 1x3 vector to define a direction vector c1-c0
%     r: the radius of the disk that is orthogonal to c1-c0, passing through c0
%     ndiv: division count to approximate a circle by a polygon, if ignored, ndiv=20
%
% output:
%     node: the 3D vertices of the disk
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

len=sqrt(sum((c0-c1).*(c0-c1)));
v0=c1-c0;

if(nargin<=2)
    r=1;
end
if(nargin<=3)
    ndiv=20;
end

dt=2*pi/ndiv;
theta=dt:dt:2*pi;
cx=r*cos(theta);
cy=r*sin(theta);
pp=[cx(:) cy(:) zeros(ndiv,1)];
node=rotatevec3d(pp,v0)+repmat(c0(:)',size(pp,1),1);
