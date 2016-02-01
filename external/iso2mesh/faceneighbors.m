function facenb=faceneighbors(t,opt)
%
% facenb=faceneighbors(t,opt)
%
% to find 4 face-neighboring elements of a tetrahedron
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%     t: tetrahedron element list, 4 columns of integers
%     opt: if opt='surface', return boundary triangle list 
%          (should be the same as the face output from v2m)
%
%          otherwise, return the element list for each element:
%          each row contains 4 numbers, representing the element
%          indices sharing triangular faces [1 2 3],[1 2 4],[1 3 4]
%          and [2 3 4] in order, where 1~4 is the node local index.
%          if the index is 0, indicating the face has no neighbor
%          (i.e. a boundary face)
%
% output:
%     facenb: see opt
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

faces=[t(:,[1,2,3]);
       t(:,[1,2,4]);
       t(:,[1,3,4]);
       t(:,[2,3,4])];
faces=sort(faces,2);
[foo,ix,jx]=unique(faces,'rows');
if(isoctavemesh)
        u=unique(jx);
        qx=u(hist(jx,u)==2);
else
        vec=histc(jx,1:max(jx));
        qx=find(vec==2);
end

nn=max(t(:));
ne=size(t,1);
facenb=zeros(size(t));

% now I need to find all repeatitive elements
% that share a face, to do this, unique('first')
% will give me the 1st element, and 'last' will
% give me the second. There will be no more than 2

% doing this is 60 times faster than doing find(jx==qx(i))
% inside a loop

if(isoctavemesh || datenum(version('-date'))>datenum('January 27 2006')) % compare to matlab 7.2
	[ujx,ii]=unique(jx,'first');
	[ujx,ii2]=unique(jx,'last');
else
	ujx=unique(jx);
	[t1,ii2]=ismember(ujx,jx);
	[t1,ii]=ismember(ujx,flipwd(jx(:)));
	ii=length(jx)-ii+1;
end

% iddup is the list of all pairs that share a common face

iddup=[ii(qx) ii2(qx)];
faceid=ceil(iddup/ne);
eid=mod(iddup,ne);
eid(eid==0)=ne;

% now rearrange this list into an element format

for i=1:length(qx)
	facenb(eid(i,1),faceid(i,1))=eid(i,2);
	facenb(eid(i,2),faceid(i,2))=eid(i,1);
end

% facenb may contain 0s, that just means the corresponding
% face is a boundary face and has no neighbor.

% if the second option is 'surface', I am going to find 
% and return surface patches only

if(nargin==2)
  if(strcmp(opt,'surface'))
	facenb=faces(find(facenb==0),:);
  else
        error(['supplied option "' opt '" is not supported.']);
  end
end
