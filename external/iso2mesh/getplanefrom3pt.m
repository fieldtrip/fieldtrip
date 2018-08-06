function [a,b,c,d]=getplanefrom3pt(plane)
%
% [a,b,c,d]=getplanefrom3pt(plane)
% 
% define a plane equation ax+by+cz+d=0 from three 3D points
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input: 
%    plane: a 3x3 matrix with each row specifying a 3D point (x,y,z)
%
% output:
%    a,b,c,d: the coefficient for plane equation ax+by+cz+d=0
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

x=plane(:,1);
y=plane(:,2);
z=plane(:,3);

% compute the plane equation a*x + b*y +c*z +d=0

a=y(1)*(z(2)-z(3))+y(2)*(z(3)-z(1))+y(3)*(z(1)-z(2));
b=z(1)*(x(2)-x(3))+z(2)*(x(3)-x(1))+z(3)*(x(1)-x(2));
c=x(1)*(y(2)-y(3))+x(2)*(y(3)-y(1))+x(3)*(y(1)-y(2));
d=-det(plane);
