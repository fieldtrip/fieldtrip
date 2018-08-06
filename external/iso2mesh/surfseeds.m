function seeds=surfseeds(node,face)
%
% seeds=surfseeds(node,face)
%
% calculate a set of interior points with each enclosed by a closed
% component of a surface
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%   node: a list of node coordinates (nn x 3)
%   face: a surface mesh triangle list (ne x 3)
%
% output:
%   seeds: the interior points coordinates for each closed-surface
%          component
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fc=finddisconnsurf(face(:,1:3));
len=length(fc);
seeds=zeros(len,3);
for i=1:len
    seeds(i,:)=surfinterior(node,fc{i});
end