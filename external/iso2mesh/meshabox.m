function [node,face,elem]=meshabox(p0,p1,opt,nodesize)
%
% [node,face,elem]=meshabox(p0,p1,opt,maxvol)
%
% create the surface and tetrahedral mesh of a box geometry
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input: 
%   p0:  coordinates (x,y,z) for one end of the box diagnoal
%   p1:  coordinates (x,y,z) for the other end of the box diagnoal
%   opt: maximum volume of the tetrahedral elements
%   nodesize: 1 or a 8x1 array, size of the element near each vertex
%
% output:
%   node: node coordinates, 3 columns for x, y and z respectively
%   face: integer array with dimensions of NB x 3, each row represents
%         a surface mesh face element 
%   elem: integer array with dimensions of NE x 4, each row represents
%         a tetrahedron 
%
% example:
%   [node,face,elem]=meshabox([2 3 2],[6 12 15],0.1,1);
%   plotmesh(node,elem,'x>4');
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
% 

if(nargin<4)
   nodesize=1;
end
[node,elem,face]=surf2mesh([],[],p0,p1,1,opt,[],[],nodesize);
elem=elem(:,1:4);
face=face(:,1:3);
