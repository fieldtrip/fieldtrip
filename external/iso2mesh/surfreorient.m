function [newnode,newface]=surfreorient(node,face)
%
% [newnode,newface]=surfreorient(node,elem)
%
% reorder nodes in a single closed surface to ensure the norms of all
% triangles are pointing outward
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2012/07/06
%
% input:
%    node: list of nodes
%    face: list of surface triangles (each row are indices of nodes of each triangle)
%
% output:
%    newnode: the output node list, in most cases it equals node
%    newface: the face list with consistent ordering
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[newnode,newface]=meshcheckrepair(node(:,1:3),face(:,1:3),'deep');
