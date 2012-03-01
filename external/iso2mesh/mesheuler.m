function [X,V,E,F]=mesheuler(face)
%
% [X,V,E,F]=mesheuler(face)
%
% Euler's charastistics of a mesh
%
% author: Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input: 
%   face: a closed surface mesh
%
% output:
%   X: Euler's charastistics
%   V: number of vertices 
%   E: number of edges
%   F: number of faces
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
% 

% mesh vertices
V=length(unique(face));

% mesh total edges
E=face(:,[1,end]);
for i=1:size(face,2)-1
	E=[E;face(:,[i,i+1])];
end
E=size(unique(sortrows(sort(E,2)),'rows'),1);

% mesh faces
F=size(face,1);

% Euler's formula, X=2-2*g and g is genus
X=V-E+F;
