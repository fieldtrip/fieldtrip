function [R,s]=rotmat2vec(u,v)
%
% [R,s]=rotmat2vec(u,v)
%
% the rotation matrix from vector u to v, satisfying R*u*s=v
%
% author: Bruno Luong
% URL:http://www.mathworks.com/matlabcentral/newsreader/view_original/827969
%
% input: 
%   u: a 3D vector in the source coordinate system;
%   v: a 3D vector in the target coordinate system;
%
% output:
%   R: a rotation matrix to transform normalized u to normalized v
%   s: a scaling factor, so that R*u(:)*s=v(:)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

s=norm(v(:))/norm(u(:));
u1=u(:)/norm(u(:));
v1=v(:)/norm(v(:));

k = cross(u1,v1);
if(~any(k)) % u and v are parallel
    R=eye(3);
    return;
end
% Rodrigues's formula:
costheta = dot(u1,v1);
R =[ 0 -k(3) k(2);
     k(3) 0 -k(1);
    -k(2) k(1) 0];
R = costheta*eye(3) + R + k*k'*(1-costheta)/sum(k.^2);
