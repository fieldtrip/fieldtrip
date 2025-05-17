function [idx, dist] = closestnode(node, p)
%
% [idx, dist]=closestnode(node,p)
%
% Find the closest point in a node list and return its index

% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    node: each row is an N-D node coordinate
%    p: a given position in the same space
%
% output:
%    idx: the index of the position in the node list that has the shortest
%         Euclidean distance to the position p
%    dist: the distances between p and each node
%
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v3 or later, see LICENSE.txt for details
%

dd = node - repmat(p, size(node, 1), 1);
[dist, idx] = min(sum(dd .* dd, 2));
