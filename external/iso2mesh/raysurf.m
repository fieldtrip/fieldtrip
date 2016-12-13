function [t,u,v,idx,xnode]=raysurf(p0,v0,node,face)
%
% [t,u,v,idx,xnode]=raysurf(p,v,node,face)
%
% perform a Havel-styled ray tracing for a triangular surface
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%   p0: list of starting points of the rays
%   v0: directional vector of the rays, 
%   node: a list of node coordinates (nn x 3)
%   face: a surface mesh triangle list (ne x 3)
%
% output:
%   t: distance from p0 to the intersection point for each surface
%      triangle, if t(i)=NaN, no intersection was found for that ray
%   u: bary-centric coordinate 1 of all intersection points
%   v: bary-centric coordinate 2 of all intersection points
%      the final bary-centric triplet is [u,v,1-u-v]
%   idx: idx lists the IDs of the face elements that intersects 
%      each ray
%   xnode: optional output, if requested, xnode gives the intersection
%      point coordinates; to compute manually, xnode=p0+repmat(t,1,3).*v0
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

len=size(p0,1);
if(len==0)
   error('p0 can not be empty');
end
if(size(node,2)<3)
   error('node must contain at least 3 columns');
end
if(size(face,2)<3)
   error('face must contain at least 3 columns');
end

if(size(v0,1)==1 | size(v0,2)==1 & len>1)
   v0=repmat(v0(:)',len,1);
end

t=zeros(len,1)*nan;
u=t;
v=t;
idx=t;

for i=1:len
   [ti,ui,vi,id]=raytrace(p0(i,:),v0(i,:),node,face);
   if(isempty(id)) continue; end
   ti=ti(id);
   tpid=find(ti>=0);
   if(isempty(tpid)) continue; end
   [tmin,tloc]=min(ti(find(ti>=0)));
   t(i)=tmin;
   u(i)=ui(id(tpid(tloc)));
   v(i)=vi(id(tpid(tloc)));
   idx(i)=id(tpid(tloc));
end

if(nargout>=5)
   xnode=p0+repmat(t,1,3).*v0;
end

