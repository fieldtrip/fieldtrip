function [node,face,elem]=meshanellip(c0,rr,tsize,maxvol)
%
% [node,face,elem]=meshanellip(c0,rr,opt)
%
% create the surface and tetrahedral mesh of an ellipsoid
%
% author: Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input: 
%   c0:  center coordinates (x0,y0,z0) of the ellipsoid
%   rr:  radii of an ellipsoid, 
%        if rr is a scalar, this is a sphere with radius rr
%        if rr is a 1x3 or 3x1 vector, it specifies the ellipsoid radii [a,b,c]
%        if rr is a 1x5 or 5x1 vector, it specifies [a,b,c,theta,phi]
%           where theta and phi are the rotation angles along z and x 
%           axes, respectively. Rotation is applied before translation.
%   tsize: maximum surface triangle size on the sphere
%   maxvol: maximu volume of the tetrahedral elements
%
% output:
%   node: node coordinates, 3 columns for x, y and z respectively
%   face: integer array with dimensions of NB x 3, each row represents
%         a surface mesh face element 
%   elem: integer array with dimensions of NE x 4, each row represents
%         a tetrahedron; if ignored, only produces the surface
%
% example:
%   [node,face,elem]=meshanellip([10,10,-10],[30,20,10,pi/4,pi/4],0.5,0.4);
%   plotmesh(node,elem,'x>10');axis equal;
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
% 

if(nargin<3)
    error('you must at least provide c0, rr and tsize, see help for details');
end

rr=rr(:)';
if(length(rr)==1)
    rr=[rr,rr,rr];
elseif(length(rr)==3)
    % do nothing
elseif(length(rr)~=5)
    error('the number of elements for rr is not correct. see help for details');
end

rmax=min(rr(1:3));
if(nargout==3)
    if(nargin==3)
	maxvol=tsize*tsize*tsize;
    end
    [node,face,elem]=meshunitsphere(tsize/rmax,maxvol/(rmax*rmax*rmax));
else
    [node,face]=meshunitsphere(tsize/rmax);
end

node=node*diag(rr(1:3));
if(length(rr)==5)
   theta=rr(4);
   phi=rr(5);
   Rz=[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];
   Rx=[1 0 0;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];
   node=(Rz*Rx*(node'))';
end
node=node+repmat(c0(:)',size(node,1),1);
