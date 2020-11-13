function [newnode,newelem]=mergemesh(node,elem,varargin)
%
% [newnode,newelem]=mergemesh(node,elem,varargin)
%
% concatenate two or more tetrahedral meshes or triangular surfaces
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input: 
%      node: node coordinates, dimension (nn,3)
%      elem: tetrahedral element or triangle surface (nn,3) to (nn,5)
%
% output:
%      newnode: the node coordinates after merging, dimension (nn,3)
%      newelem: tetrahedral element or surfaces after merging (nn,4) or (nhn,5)
%
% note: you can call meshcheckrepair for the output newnode and
% newelem to remove the duplicated nodes or elements. mergemesh does
% detect self-intersecting elements when merging; to remove self-intersecting
% elements, you need to use mergesurf().
%
% example:
%
%   [node1,face1,elem1]=meshabox([0 0 0],[10 10 10],1,1);
%   [node2,face2,elem2]=meshasphere([5 5 13.1],3,0.3,3);
%   [newnode,newelem]=mergemesh(node1,elem1,node2,elem2);
%   plotmesh(newnode,newelem);
%   figure;
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

X=mesheuler(newelem);

if(size(newelem,2)==4)
   if(X>=0)
      newelem(:,end+1)=1;
   end
end
if(size(newelem,2)==3)
   newelem(:,end+1)=1;
end
for i=1:2:len
   no=varargin{i};
   el=varargin{i+1};
   baseno=size(newnode,1);
   if(size(no,2)~=size(newnode,2))
        error('input node arrays have inconsistent columns');
   end
   if(size(el,2)==5 | size(el,2)==4)
        el(:,1:4)=el(:,1:4)+baseno;
	if(size(el,2)==4 & X>=0)
	   el(:,5)=1+(i+1)/2;
	end
   	newnode=[newnode;no];
	newelem=[newelem;el];
   elseif(size(el,2)==3 & size(newelem,2)==4)
        el(:,1:3)=el(:,1:3)+baseno;
	if(size(el,2)==3)
	   el(:,4)=1+(i+1)/2;
	end
   	newnode=[newnode;no];
	newelem=[newelem;el];
   else
        error('input element arrays have inconsistent columns');
   end
end
