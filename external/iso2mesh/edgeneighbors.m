function edgenb=edgeneighbors(t,opt)
%
% edgenb=edgeneighbors(t,opt)
%
% to find neighboring triangular elements in a triangule surface
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%     t: a triangular surface element list, 3 columns of integers
%     opt: if opt='general', return the edge neighbors for a general
%          triangular surface: each edge can be shared by more than 2
%          triangles; if ignored, we assume all triangles are shared by no
%          more than 2 triangles.
%
% output:
%     edgenb: if opt is not supplied, edgenb is a size(t,1) by 3 array with
%     each element being the triangle ID of the edge neighbor of that
%     triangle. For each row, the order of the neighbors is listed as those
%     sharing edges [1 2], [2 3] and [3 1] between the triangle nodes.
%
%     when opt='general', edgenb is a cell array with a length of size(t).
%     each member of the cell array is a list of edge neighbors (the order 
%     is not defined).
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

edges=[t(:,[1,2]);
       t(:,[2,3]);
       t(:,[3 1])];
edges=sort(edges,2);
[foo,ix,jx]=unique(edges,'rows');

if(nargin==2)
  if(strcmp(opt,'general'))
        ne=size(t,1);
        edgenb=cell(ne,1);
        for i=1:ne
            % this is very slow, need to be optimized
            nb=unique(mod([find(jx==jx(i) | jx==jx(i+ne) | jx==jx(i+2*ne))]',ne),'first');
            nb(nb==0)=ne;
            edgenb{i}=nb(nb~=i);
        end
        return;
  else
        error(['supplied option "' opt '" is not supported.']);
  end
end

if(isoctavemesh)
        u=unique(jx);
        qx=u(hist(jx,u)==2);
else
        vec=histc(jx,1:max(jx));
        qx=find(vec==2);
end

nn=max(t(:));
ne=size(t,1);
edgenb=zeros(size(t));

% now I need to find all repeatitive elements
% that share a face, to do this, unique('first')
% will give me the 1st element, and 'last' will
% give me the second. There will be no more than 2

% doing this is 60 times faster than doing find(jx==qx(i))
% inside a loop

[ujx,ii]=unique(jx,'first');
[ujx,ii2]=unique(jx,'last');

% iddup is the list of all pairs that share a common face

iddup=[ii(qx) ii2(qx)];
faceid=ceil(iddup/ne);
eid=mod(iddup,ne);
eid(eid==0)=ne;

% now rearrange this list into an element format

for i=1:length(qx)
	edgenb(eid(i,1),faceid(i,1))=eid(i,2);
	edgenb(eid(i,2),faceid(i,2))=eid(i,1);
end

% edgenb may contain 0s, that just means the corresponding
% face is a boundary face and has no neighbor.

% if the second option is 'surface', I am going to find 
% and return surface patches only


