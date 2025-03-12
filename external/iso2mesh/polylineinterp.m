function [idx, weight, newnodes] = polylineinterp(polylen, len, nodes)
%
% [idx, weight]=polylineinterp(polylen, len)
% [idx, weight, newnodes]=polylineinterp(polylen, len, nodes)
%
% Find the polyline segment indices and interpolation weights for a
% specified total length or a set of lengths
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    polylen: a 1D vector sequentially recording the length of each segment
%         of a polyline, the first number is the length of the 1st segment,
%         and so on
%    len: a single scalar, or a vector of scalars, specifying the total
%         length
%    nodes: if nodes is an array with a row-number equal to length(polylen)+1,
%         we assume each row defines a coordinate for the nodes along the
%         polyline
%
% output:
%    idx: the indices of the polyline segments, starting from 1, where each
%         length defined in len ends; if len> sum(polylen), nan is
%         returned; if len<0, the weight will be a negative value.
%    weight: the interpolation weight between 0-1 towards the end node
%         of the containing segment; the weight for the start-node is 1-weight
%    newnodes: the interpolated node positions at the end of the len
%
% example:
%    lineseg=[2,2,1,7,10];
%    [idx, weight]=polylineinterp(lineseg, [3, 12, 7])
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v3 or later, see LICENSE.txt for details
%

cumlen = [0 cumsum(polylen(:)')];
idx = nan * ones(size(len));
weight = zeros(size(len));

if (nargin >= 3 && nargout >= 3)
    if (size(nodes, 1) == 1)
        nodes = nodes.';
    end
    if (size(nodes, 1) ~= length(polylen) + 1)
        error('the row number of the nodes input must be 1 more than the length of polylen');
    end
    newnodes = zeros(length(len), size(nodes, 2));
end

for i = 1:length(len)
    pos = histc(len(i), cumlen);
    if (any(pos == 1))
        idx(i) = find(pos);
        if (idx(i) == length(cumlen))
            idx(i) = idx(i) - 1;
            weight(i) = 1;
            newnodes(i, :) = nodes(end, :);
        elseif (idx(i) <= length(polylen))
            weight(i) = (len(i) - cumlen(idx(i))) / polylen(idx(i));
            if (nargin >= 3 && nargout >= 3)
                newnodes(i, :) = (1 - weight(i)) * nodes(idx(i), :) + weight(i) * nodes(idx(i) + 1, :);
            end
        end
    end
end

idx(idx > length(polylen)) = nan;
