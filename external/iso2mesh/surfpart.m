function elist=surfpart(f,loopedge)
%
% elist=surfpart(f,loopedge)
%
% partition a triangular surface using a closed loop defined by existing edges
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2012/02/09
%
% input:
%      f: input, surface face element list, dimension (be,3)
%      loopedge: a 2-column array, specifying a closed loop in CCW order
%
% output:
%      elist: list of triangles that is enclosed by the loop
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

elist=[];
if(isempty(f) || isempty(loopedge))
    return;
end

if(size(f,2)==3)
    edges=[f(:,[1,2]);
           f(:,[2,3]);
           f(:,[3,1])];             % create all the edges
elseif(size(f,2)==4)
    edges=[f(:,[1,2]);
           f(:,[2,3]);
           f(:,[3,4]);
           f(:,[4,1])];             % create all the edges
else
    error('surfpart only supports triangular and quadrilateral elements');
end

[elist,front]=advancefront(edges,loopedge);
while(~isempty(front))
	[elist0,front0]=advancefront(edges,front);
	elist=unique([elist;elist0]);
	front=front0;
end
