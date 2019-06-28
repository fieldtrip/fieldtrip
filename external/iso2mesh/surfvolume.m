function vol=surfvolume(node,face,option)
%
% vol=surfvolume(node,face,option)
%
% calculate the enclosed volume for a closed surface
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    node:  node coordinates
%    face:  surface triangle list
%
% output:
%    vol:   total volume of the enclosed space
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

face=face(:,1:3);

ed=surfedge(face);
if(~isempty(ed))
   error('open surface is detected, you have to close it first, consider meshcheckrepair() with meshfix option');
end

[no,el]=fillsurf(node,face);

vol=elemvolume(no,el);
vol=sum(vol);
