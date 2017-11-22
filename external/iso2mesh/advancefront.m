function [elist,nextfront]=advancefront(edges,loop,elen)
%
% [elist,nextfront]=advancefront(edges,loop,elen)
%
% advance an edge-front on an oriented surface to the next separated by
% one-element width
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2012/02/09
%
% input:
%      edges: edge list of an oriented surface mesh, must be in CCW order
%      loop: a 2-column array, specifying a closed loop in CCW order
%      elen: node number inside each element, if ignored, elen is set to 3
%
% output:
%      elist: list of triangles that is enclosed between the two
%             edge-fronts
%      nextfront: a new edge loop list representing the next edge-front
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

elist=[];
nextfront=[];

if(nargin<3) elen=3; end

[hasedge, loc]=ismember(loop,edges,'rows');

if(~all(hasedge))
   error('loop edge is not defined in the mesh');
end

nodenum=size(edges,1)/elen;

elist=unique(mod(loc-1,nodenum))+1;
nextfront=edges(elist,:);
for i=1:elen-1
    nextfront=[nextfront;edges(elist+nodenum*i,:)];
end
nextfront=setdiff(nextfront,loop,'rows');

% remove reversed edge pairs
[flag,loc]=ismember(nextfront,nextfront(:,[2 1]),'rows');
id=find(flag);
if(~isempty(id))
    delmark=flag;
    delmark(loc(find(loc>0)))=1;
    nextfront(find(delmark),:)=[];
end
nextfront=nextfront(:,[2 1]); % reverse this loop, as it is reversed to the input loop
