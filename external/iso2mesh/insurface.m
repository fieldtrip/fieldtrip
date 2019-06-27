function tf=insurface(node,face,points)
%
% tf=innersurf(node,face,points)
%
% test if a set of 3D points is located inside a 3D triangular surface
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    node:  node coordinates
%    face:  surface triangle list
%    points: a set of 3D points (Nx3 array)
%
% output:
%    tf: a vector with the same length of points, 
%        a value of 1 means the point is inside of the surface, and
%        a value of 0 means the point is outside of the surface.
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[no,el]=fillsurf(node,face);
tf=tsearchn(no,el,points);

tf(~isnan(tf))=1;
tf(isnan(tf)) =0;
