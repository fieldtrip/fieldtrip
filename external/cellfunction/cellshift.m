function [x] = cellshift(x, shift, dim, maxshift)

% CELLSHIFT(X, SHIFT, DIM) 

if numel(shift)>1
  x = x(:)';     
  y = cell(numel(shift),numel(x));
  for k = 1:numel(shift)
    %y(k,:) = cellshift(x, shift(k), dim, max(abs(shift)));
    y(k,:) = cellshift(x, shift(k), dim, [abs(min(shift)) abs(max(shift))]);
  end
  for k = 1:size(y,2)
    y{1,k} = cat(1,y{:,k});
  end
  x = y(1,:);
  return;
end

if nargin<4
  maxshift = max(abs(shift));
else
  if any(abs(shift))>abs(maxshift)
    error('the value for maxshift should be >= shift');
  end
end


if nargin<3
  dim = find(size(x{1})>1, 1, 'first');
end

maxshift = abs(maxshift);
if numel(maxshift)==1, maxshift = maxshift([1 1]); end 

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1)
  error('incorrect input for cellshift');
end

n    = numel(x);
nsmp = cellfun('size', x, dim);
beg1 = ones(1,n) + shift + maxshift(1);
end1 = nsmp      + shift - maxshift(2);

switch dim
case 1
  for k = 1:n
    x{k} = x{k}(beg1(k):end1(k),:);
  end
case 2
  for k = 1:n
    x{k} = x{k}(:,beg1(k):end1(k));
  end
otherwise
  error('dimensionality of >2 is not supported');
end
