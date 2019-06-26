function [no,el]=fillsurf(node,face)
%
% [no,el]=fillsurf(node,face)
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

ISO2MESH_TETGENOPT='-YY';
[no,el]=surf2mesh(node,face,[],[],1,1);
