function [newnode,newelem]=mergesurf(node,elem,varargin)
%
% [newnode,newelem]=mergesurf(node1,elem1,node2,elem2,...)
%
% merge two or more triangular meshes and split intersecting elements
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input:
%      node: node coordinates, dimension (nn,3)
%      elem: tetrahedral element or triangle surface (nn,3)
%
% output:
%      newnode: the node coordinates after merging, dimension (nn,3)
%      newelem: tetrahedral element or surfaces after merging (nn,4) or (nhn,5)
%
% note: you can call meshcheckrepair for the output newnode and
% newelem to remove the duplicated nodes or elements
%
% example:
%
%   [node1,face1,elem1]=meshabox([0 0 0],[10 10 10],1,1);
%   [node2,face2,elem2]=meshasphere([5 5 10],3,0.3,3);
%   [newnode,newface]=mergemesh(node1,face1,node2,face2);
%   plotmesh(newnode,newface,'x>5');
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

len=length(varargin);
newnode=node;
newelem=elem;
if(len>0 & mod(len,2)~=0)
   error('you must give node and element in pairs');
end
for i=1:2:len
   no=varargin{i};
   el=varargin{i+1};
   [newnode,newelem]=surfboolean(newnode,newelem,'all',no,el);
end
