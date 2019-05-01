function [inface, outface]=innersurf(node,face,outface)
%
% outface=innersurf(node,face,outface)
%
% extract the interior triangles (shared by two enclosed compartments) of a complex surface
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    node:  node coordinates
%    face:  surface triangle list
%    outface: (optional) the exterior triangle list, if not given, will
%           be computed using outersurf().
%
% output:
%    inface: the collection of interior triangles of the surface mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin<3)
  outface=outersurf(node,face);
end

[I,J]=ismember(sort(face,2),sort(outface,2),'rows');

inface=face(I==0,:);
