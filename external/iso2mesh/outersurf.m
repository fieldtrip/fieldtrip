function outface=outersurf(node,face)
%
% outface=outersurf(node,face)
%
% extract the out-most shell of a complex surface mesh
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    node:  node coordinates
%    face:  surface triangle list
%
% output:
%    outface: the out-most shell of the surface mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

face=face(:,1:3);

ed=surfedge(face);
if(~isempty(ed))
   error('open surface is detected, you have to close it first, consider meshcheckrepair() with meshfix option');
end

[no,el]=fillsurf(node,face);

outface=volface(el);

[no,outface]=removeisolatednode(no,outface);
maxfacenode=max(outface(:));

[I,J]=ismember(round(no(1:maxfacenode,:)*1e10),round(node*1e10),'rows');

% if(sum(I)~=maxfacenode)
%     error('mesh tessellation failed');
% end
outface=J(outface);
[ii,jj]=find(outface==0);
outface(ii,:)=[];
