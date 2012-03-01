function newedge=orderloopedge(edge)
%
% [newedge]=orderloopedge(edge)
%
% order the node list of a simple loop based on connection sequence
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
% date:   2008/05
%
% input: 
%        edge: a loop consisted by a sequence of edges, each row 
%              is an edge with two integers: start/end node index
%
% output:
%        newedge: reordered edge node list
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

% this subroutine can not process bifercation

ne=size(edge,1);
newedge=zeros(size(edge));
newedge(1,:)=edge(1,:);

for i=2:ne
  [row,col]=find(edge(i:end,:)==newedge(i-1,2));
  if(length(row)==1) % loop node
     newedge(i,:)=[newedge(i-1,2),edge(row+i-1,3-col)];
     edge([i,row+i-1],:)=edge([row+i-1,i],:);
  elseif (length(row)>=2)
     error('bifercation is found,exit');
  elseif (length(row)==0)
     error(['open curve at ' num2str(edge(i-1,2)) ]);
  end
end
     
