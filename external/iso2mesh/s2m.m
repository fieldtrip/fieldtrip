function [node,elem,face]=s2m(v,f,keepratio,maxvol)
%
% [node,elem,face]=s2m(v,f,keepratio,maxvol)
%
% volumetric mesh generation from a closed surface, shortcut for surf2mesh
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% inputs and outputs are similar to those defined in surf2mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

p0=min(v(:,1:3));
p1=max(v(:,1:3));
[node,elem,face]=surf2mesh(v,f,p0,p1,keepratio,maxvol,[],[]);

