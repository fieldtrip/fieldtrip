function newpt=rotatevec3d(pt,v1,u1,p0)
%
% newpt=rotatevec3d(pt,v1,u1,p0)
%
% rotate 3D points from one Cartesian coordinate system to another
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input: 
%   pt: 3D points defined in a standard Cartesian system where a unitary
%       z-vector is (0,0,1), 3 columns for x, y and z 
%   v1: the unitary z-vector for the target coordinate system
%   u1: the unitary z-vector for the source coordinate system, if ignored,
%       u1=(0,0,1) 
%   p0: offset of the new coordinate system, if ignored, p0=(0,0,0)
%
% output:
%   newpt: the transformed 3D points
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin<2)
    error('you must give at least pt and v1');
end
if(nargin==2)
    u1=[0,0,1];
end
if(nargin<=3)
    p0=[0,0,0];
end

u1=u1/norm(u1);
v1=v1/norm(v1);

[R,s]=rotmat2vec(u1,v1);
newpt=(R*pt'*s)';

if(nargin>3)
  p0=p0(:)';
  newpt=newpt+repmat(p0,size(newpt,1),1);
end