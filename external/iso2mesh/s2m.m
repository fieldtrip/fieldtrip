function [node,elem,face]=s2m(v,f,keepratio,maxvol,method)
%
% [node,elem,face]=s2m(v,f,keepratio,maxvol)
%
% volumetric mesh generation from a closed surface, shortcut for surf2mesh
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% inputs and outputs are similar to those defined in surf2mesh
%
% if method='cgalpoly', s2m will call cgals2m and keepratio should be a 
% structure (as the 'opt' input in cgals2m)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

p0=min(v(:,1:3));
p1=max(v(:,1:3));
if(nargin>=5)
  if(strcmp(method,'cgalpoly'))
    [node,elem,face]=cgals2m(v,f,keepratio,maxvol);
    return;
  end
end
[node,elem,face]=surf2mesh(v,f,p0,p1,keepratio,maxvol,[],[]);

