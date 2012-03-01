function plane=surfplane(node,face)
%
% plane=surfplane(node,face)
%
% plane equation coefficients for each face in a surface
%
% author: Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input:
%   node: a list of node coordinates (nn x 3)
%   face: a surface mesh triangle list (ne x 3)
%
% output:
%   plane: a (ne x 4) array, in each row, it has [a b c d]
%        to denote the plane equation as "a*x+b*y+c*z+d=0"
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

AB=node(face(:,2),1:3)-node(face(:,1),1:3);
AC=node(face(:,3),1:3)-node(face(:,1),1:3);

N=cross(AB',AC')';
d=-dot(N',node(face(:,1),1:3)')';
plane=[N,d];
