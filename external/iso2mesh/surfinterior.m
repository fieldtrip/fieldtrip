function [pt,p0,v0,t,idx]=surfinterior(node,face)
%
% [pt,p0,v0,t,idx]=surfinterior(node,face)
%
% identify a point that is enclosed by the (closed) surface
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%   node: a list of node coordinates (nn x 3)
%   face: a surface mesh triangle list (ne x 3)
%
% output:
%   pt: the interior point coordinates [x y z]
%   p0: ray origin used to determine the interior point
%   v0: the vector used to determine the interior point
%   t : ray-tracing intersection distances (with signs) from p0. the
%       intersection coordinates can be expressed as p0+t(i)*v0
%   idx: index to the face elements that intersect with the ray, order
%       match that of t
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

pt=[];
len=size(face,1);
for i=1:len
   p0=mean(node(face(i,1:3),:));
   plane=surfplane(node,face(i,:));
   v0=plane(1:3);

   [t,u,v]=raytrace(p0,v0,node,face(:,1:3));

   idx=find(u>=0 & v>=0 & u+v<=1.0 & ~isinf(t));
   [ts, uidx]=unique(sort(t(idx)));
   if(~isempty(ts) && mod(length(ts),2)==0)
       ts=reshape(ts,[2 length(ts)/2]);
       tdiff=ts(2,:)-ts(1,:);
       [maxv,maxi]=max(tdiff);
       pt=p0+v0*(ts(1,maxi)+ts(2,maxi))*0.5;
       idx=idx(uidx);
       t=t(idx);
       break;
   end
end
