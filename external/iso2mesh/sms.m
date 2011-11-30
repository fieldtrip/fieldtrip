function newnode=sms(node,face,iter,alpha,method)
%
% newnode=sms(node,face,iter,useralpha,method)
%
% simplified version of surface mesh smoothing
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2009/10/21
%
% input:
%    node:  node coordinates of a surface mesh
%    face:  face element list of the surface mesh
%    iter:  smoothing iteration number
%    alpha: scaler, smoothing parameter, v(k+1)=alpha*v(k)+(1-alpha)*mean(neighbors)
%    method: same as in smoothsurf, default is 'laplacianhc'
%
% output:
%    newnode: output, the smoothed node coordinates
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin<5)
   method='laplacianhc';
end
if(nargin<4)
   if(nargin<3)
      iter=10;
   end
   alpha=0.5;
end

conn=meshconn(face,size(node,1));
newnode=smoothsurf(node(:,1:3),[],conn,iter,alpha,method,alpha);
