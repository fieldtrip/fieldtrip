function [t,u,v,idx]=raytrace(p0,v0,node,face)
%
% [t,u,v,idx]=raytrace(p0,v0,node,face)
%
% perform a Havel-styled ray tracing for a triangular surface
%
% author: Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input:
%   p0: starting point coordinate of the ray
%   v0: directional vector of the ray
%   node: a list of node coordinates (nn x 3)
%   face: a surface mesh triangle list (ne x 3)
%
% output:
%   t: signed distance from p to the intersection point for each surface
%      triangle, if ray is parallel to the triangle, t is set to Inf
%   u: bary-centric coordinate 1 of all intersection points
%   v: bary-centric coordinate 2 of all intersection points
%      the final bary-centric triplet is [u,v,1-u-v]
%   idx: optional output, if requested, idx lists the IDs of the face
%      elements that intersects the ray; users can manually calc idx by
%
%      idx=find(u>=0 & v>=0 & u+v<=1.0 & ~isinf(t));
%
% Reference: 
%  [1] J. Havel and A. Herout, "Yet faster ray-triangle intersection (using 
%          SSE4)," IEEE Trans. on Visualization and Computer Graphics,
%          16(3):434-438 (2010)
%  [2] Q. Fang, "Comment on 'A study on tetrahedron-based inhomogeneous 
%          Monte-Carlo optical simulation'," Biomed. Opt. Express, (in
%          press)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

p0=p0(:)';
v0=v0(:)';

AB=node(face(:,2),1:3)-node(face(:,1),1:3);
AC=node(face(:,3),1:3)-node(face(:,1),1:3);

N=cross(AB',AC')';
d=-dot(N',node(face(:,1),1:3)')';

Rn2=1./sum((N.*N)')';

N1=cross(AC',N')'.*repmat(Rn2,1,3);
d1=-dot(N1',node(face(:,1),1:3)')';

N2=cross(N',AB')'.*repmat(Rn2,1,3);
d2=-dot(N2',node(face(:,1),1:3)')';

den=(v0*N')';
t=-(d+(p0*N')');
P=(p0'*den'+v0'*t')';
u=dot(P',N1')'+den.*d1;
v=dot(P',N2')'+den.*d2;

idx=find(den);
den(idx)=1./den(idx);

t=t.*den;
u=u.*den;
v=v.*den;

 % if den==0, ray is parallel to triangle, set t to infinity
t(find(den==0))=Inf;

if(nargout>=4)
    idx=find(u>=0 & v>=0 & u+v<=1.0 & ~isinf(t));
end
