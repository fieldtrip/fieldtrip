function islands = bwislands(img)
%
% islands=bwislands(img)
%
% return the indices of non-zero elements in a 2D or 3D image
% grouped by connected regions in a cell array
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    img: a 2D or 3D array
% output:
%    islands: a cell array, each cell records the indices
%         of the non-zero elements in img for a connected
%         region (or an island)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

img = logical(1 - img);
idx = find(1 - img(:));
islands = {};

count = 1;
out = cell(1, ndims(img));
while (~isempty(idx))
    [out{:}] = ind2sub(size(img), idx(1));
    imgnew = imfill(img, cell2mat(out));
    islands{count} = find(imgnew ~= img);
    count = count + 1;
    img = imgnew;
    idx = find(1 - img(:));
end
