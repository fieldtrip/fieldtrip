function [f maxsize]=maxsurf(facecell,node)
%
% [f maxsize]=maxsurf(facecell,node)
%
% return the surface with the maximum element number or  
% total area from a cell arry of surfaces
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    facecell: a cell array, each element is a face array
%    node: optional, node list, if given, the output is the
%          surface with the largest surface area.
%
% output:
%    f: the surface data (node indices) for the surface with the 
%       most elements (or largest area when node is given)
%    maxsize: if node is not given, maxisize is row number of f;
%             otherwise, maxsize is the total area of f
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

maxsize=-1;
maxid=-1;

if(nargin==2)
    areas=zeros(1,length(facecell));
    for i=1:length(facecell)
       areas(i)=sum(elemvolume(node(:,1:3),facecell{i}));
    end
    [maxsize,maxid]=max(areas);
    f=facecell{maxid};
    return;
else
    for i=1:length(facecell)
	if(length(facecell{i})>maxsize)
		maxsize=length(facecell{i});
		maxid=i;
	end
    end
    f=[];
    if(maxid>0)
	f=facecell{maxid};
    end
end
