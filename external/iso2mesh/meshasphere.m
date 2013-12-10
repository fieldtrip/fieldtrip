function [node,face,elem]=meshasphere(c0,r,tsize,maxvol)
%
% [node,face,elem]=meshasphere(c0,r,tsize,maxvol)
%
% create the surface and tetrahedral mesh of a sphere
%
% author: Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input: 
%   c0:  center coordinates (x0,y0,z0) of the sphere
%   r:   radius of the sphere
%   tsize: maximum surface triangle size on the sphere
%   maxvol: maximu volume of the tetrahedral elements
%
% output:
%   node: node coordinates, 3 columns for x, y and z respectively
%   face: integer array with dimensions of NB x 3, each row represents
%         a surface mesh face element 
%   elem: integer array with dimensions of NE x 4, each row represents
%         a tetrahedron 
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
% 

if(nargin<3)
    error('you must at least provide c0, r and tsize');
end
if(nargin==3)
    maxvol=tsize*tsize*tsize;
end
if(nargout==3)
    [node,face,elem]=meshunitsphere(tsize/r,maxvol/(r*r*r));
else
    [node,face]=meshunitsphere(tsize/r);
end

node=node*r+repmat(c0(:)',size(node,1),1);
