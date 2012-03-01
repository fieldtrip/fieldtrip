function f=maxsurf(facecell)
%
% f=maxsurf(facecell)
%
% return the surface with the maximum element number (not 
% necessarily in area) from a cell arry of surfaces
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%    facecell: a cell array, each element is a face array
%
% output:
%    f: the surface data (node indices) for the surface with most elements
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

maxsize=-1;
maxid=-1;

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
