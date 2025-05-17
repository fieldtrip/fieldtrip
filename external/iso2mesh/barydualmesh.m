function [newnode, newelem] = barydualmesh(node, elem, flag)
%
% [newnode,newelem]=barydualmesh(node,elem)
%
% generate barycentric dual-mesh by connecting edge, face and elem centers
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    node: list of input mesh nodes
%    elem: list of input mesh elements (each row are indices of nodes of each element)
%    flag: if is 'cell', output newelem as cell arrays (each has 1x4 nodes)
%
% output:
%    newnode: all new nodes in the barycentric dual-mesh (made of edge/face/elem centers)
%    newelem: the indices of the face nodes for each original tet element
%
% example:
%    [node,elem]=meshgrid6([0 60],[0 60],[0 60]);
%    [newnode,newelem]=barydualmesh(node,elem,'cell');
%    plotmesh(newnode,newelem);
%    hold on; plotmesh(node,[],elem,'facecolor','none','edgecolor','b')
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[enodes, eidx] = highordertet(node, elem);               % compute edge-centers

[fnodes, fidx] = elemfacecenter(node, elem);             % compute face-centers
c0 = meshcentroid(node, elem(:, 1:min(size(elem, 2), 4))); % compute elem-centers

% concatenated new nodes and their indices
newnode = [enodes; fnodes; c0];
newidx = [eidx, fidx + size(enodes, 1), (1:size(elem, 1))' + (size(enodes, 1) + size(fnodes, 1))];

newelem = [
           1 8 11 7
           2 7 11 9
           3 9 11 8
           4 7 11 10
           5 8 11 10
           6 9 11 10
          ];
newelem = newelem';
newelem = newidx(:, newelem(:));
newelem = reshape(newelem, [size(elem, 1) 4 6]);
newelem = permute(newelem, [1 3 2]);
newelem = reshape(newelem, [size(elem, 1) * 6 4]);
if (nargin > 2 && ischar(flag) && strcmp(flag, 'cell'))
    newelem = num2cell(newelem, 2);
end
