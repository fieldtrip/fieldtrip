function [len, node, inputreversed] = polylinelen(node, p0, p1, pmid)
%
% [len, node, inputreversed]=polylinelen(node, p0, p1, pmid)
%
% Calculate the polyline line segment length vector in sequential order
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    node: an N x 3 array defining each vertex of the polyline in
%          sequential order
%    p0:(optional) a given node to define the start of the polyline, if not
%         defined, start position is assumed to be 1st node
%    p1:(optional) a given node to define the end of the polyline, if not
%         defined, end position is assumed to be last node
%    pmid:(optional) a given node sits between p0 and p1, if not
%         defined, index of the middle (floored) node is used
%
% output:
%    len: the length of each segment between the start and the end points
%    node: the node list between the start and end points of the polyline
%    inputreversed: if 1, the input node is reversed from p0 to pmid to p1
%
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v3 or later, see LICENSE.txt for details
%

if (nargin < 3)
    p1 = size(node, 1);
    if (nargin < 2)
        p0 = 1;
    end
end

if (nargin < 4)
    pmid = floor((p0 + p1) * 0.5);
end

if (size(p0, 2) == 3)
    p0 = closestnode(node, p0);
end
if (size(p1, 2) == 3)
    p1 = closestnode(node, p1);
end
if (size(pmid, 2) == 3)
    pmid = closestnode(node, pmid);
end

if (p0 < pmid && pmid < p1)
    inputreversed = 0;
    node = node(p0:p1, :);
elseif (p0 < pmid && p1 < pmid)
    inputreversed = (min(p0, p1) == p0);
    node = node([min(p0, p1):-1:1 end:-1:max(p0, p1)], :);
    if (~inputreversed)
        node = flipud(node);
    end
elseif (p0 > pmid && pmid > p1)
    inputreversed = 1;
    node = node(p0:-1:p1, :);
elseif (p0 > pmid && p1 > pmid)
    inputreversed = (max(p0, p1) == p1);
    node = node([max(p0, p1):end 1:min(p0, p1)], :);
    if (inputreversed)
        node = flipud(node);
    end
end

len = node(1:end - 1, :) - node(2:end, :);
len = sqrt(sum(len .* len, 2));
