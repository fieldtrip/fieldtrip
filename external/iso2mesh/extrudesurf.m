function [node,face]=extrudesurf(no,fc,vec)
%
% [node,face]=extrudesurf(no,fc,vec)
% 
% create a enclosed surface mesh by extruding an open surface
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%
% output:
%      node: 3D node coordinates for the generated surface mesh
%      face: triangular face patches of the generated surface mesh, each
%           row represents a triangle denoted by the indices of the 3 nodes
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

nlen=size(no,1);
if(length(vec)>1)
   node=[no; no+repmat(vec(:)', nlen,1)];
else
   node=[no; no+vec*nodesurfnorm(no, fc)];
end

face=[fc; fc+nlen];

edge=surfedge(fc);
sideface=[edge edge(:,1)+nlen; edge+nlen edge(:,2)];
face=[face; sideface];


