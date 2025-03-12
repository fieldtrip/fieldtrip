function node = orthdisk(c0, c1, r, ndiv, v1, angle0)
%
% node=orthdisk(c0,c1,r,ndiv)
%
% Defining a 3D disk that is orthogonal to the vector c1-c0
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%     c0: a 1x3 vector for the origin
%     c1: a 1x3 vector to define a direction vector c1-c0
%     r: the radius of the disk that is orthogonal to c1-c0, passing through c0
%     ndiv: division count to approximate a circle by a polygon, if ignored, ndiv=20
%     v1: a 1x3 vector specifying the first point on the output 3D disk. if
%         v1 is not perpendicular to c1-c0, the disk rotation axis is
%         changed to cross(v1,cross(c1-c0,v1));
%     angle0: angle0 represents the angle (in radian) of the 1st point in
%         the 3D disk if placed on the x-y plane.
%
% output:
%     node: the 3D vertices of the disk
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

len = sqrt(sum((c0 - c1) .* (c0 - c1)));
v0 = c1 - c0;

if (nargin >= 5)
    vt = cross(v0, v1);
    if (abs(dot(v0(:), v1(:))) > 1e-5) % input is not orthogonal
        v0 = cross(v1, vt);
    end
end

if (nargin < 6)
    angle0 = 0;
end

if (nargin <= 2)
    r = 1;
end
if (nargin <= 3)
    ndiv = 20;
end

dt = 2 * pi / ndiv;
theta = angle0 + dt:dt:2 * pi + angle0;
cx = r * cos(theta);
cy = r * sin(theta);
pp = [cx(:) cy(:) zeros(ndiv, 1)];
node = rotatevec3d(pp, v0) + repmat(c0(:)', size(pp, 1), 1);
