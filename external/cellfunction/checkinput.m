function [catdim] = checkinput(x)

if ~iscell(x)
  error('input should be a cell-array');
end

siz    = size2(x, [], 'cell');
ndim   = size(siz, 2);

if ndim>2
  error('number of dimensions within a cell > 2 is currently not supported');
end

if length(x)<numel(x)
  error('currently only a 1xN or Nx1 cell-array is supported in the input');
end

same   = false(1,ndim);
for k = 1:ndim
  same(k) = all(siz(1,k)==siz(:,k));
end

if sum(same)==1
  % data need to be concatenated across the other dimension
  catdim = ~same;
elseif sum(same)==0
  % data cannot be concatenated across any of the dimensions
  error('the cells in the input array cannot be concatenated across any of the dimensions in the cells');
elseif sum(same)>1
  % more than one possibility 
  catdim = same;
end

