function [y] = cellshift(x, shift, dim, maxshift, flag)

% CELLSHIFT(X, SHIFT, DIM) 

if nargin<5 || isempty(flag)
  flag = 'overlap'; % can also be 'zero' or 'nan'
  % overlap: only return the data samples that fully overlap for all shifts
  % zero: return 0's for the samples at which no data is present
  % nan: return nans for the samples at which no data is present
end

if nargin<4 || isempty(maxshift)
  maxshift = max(abs(shift));
else
  if any(abs(shift(:)))>abs(maxshift)
    error('the value for maxshift should be >= shift');
  end
end

if nargin<3 || isempty(dim)
  dim = find(size(x{1})>1, 1, 'first');
end
  
siz = size(shift);
if siz(1)==1 && siz(2)==1
  % this is OK, continue below
elseif siz(1)==1 && siz(2)>1
  % call the function recursively, and concatenate, same set of shifts for
  % all channels
  x = x(:)';     
  y = cell(numel(shift),numel(x));
  for k = 1:numel(shift)
    %y(k,:) = cellshift(x, shift(k), dim, max(abs(shift)));
    y(k,:) = cellshift(x, shift(k), dim, [abs(min(shift)) abs(max(shift))], flag);
  end
  for k = 1:size(y,2)
    y{1,k} = cat(1,y{:,k});
  end
  y = y(1,:);
  return;
elseif siz(1)>1 && siz(2)>=1
  % different shifts, call the function recursively, and concatenate.
  % note that the order of the rows in the output is different from when
  % the same shift applies to all channels
  edges =[abs(nanmin(shift(:))) abs(nanmax(shift(:)))];
  y     = cell(siz(1),numel(x));
  for k = 1:siz(1)
    tmpshift = shift(k,:);
    tmpshift(~isfinite(tmpshift)) = [];
    if isempty(tmpshift)
      tmpshift = 0;
    end
    y(k,:) = cellshift(cellrowselect(x, k), tmpshift, dim, edges, flag);
  end
  for k = 1:size(y,2)
    y{1,k} = cat(1,y{:,k});
  end
  y = y(1,:);
  return;
end

%-----------------------------------------------------
maxshift = abs(maxshift);
if numel(maxshift)==1, maxshift = maxshift([1 1]); end 

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1)
  error('incorrect input for cellshift');
end

n    = numel(x);
nsmp = cellfun('size', x, dim);

switch flag
  case 'overlap'
    beg1 = ones(1,n) + shift + maxshift(1);
    end1 = nsmp      + shift - maxshift(2);

    beg2 = ones(size(beg1));
    end2 = end1-beg1+1;
    
    nsmp = nsmp - sum(maxshift);
    
  case {'zero' 'nan'}
    beg1 = ones(1,n);
    end1 = nsmp;
    
    beg2 = ones(1,n) - shift + maxshift(2);
    end2 = nsmp      - shift + maxshift(2); 
    
    nsmp = nsmp + sum(maxshift);
end

usenans = strcmp(flag, 'nan');

switch dim
case 1
  y   = cell(size(x));
  for k = 1:n
    tmp = zeros(nsmp(k), size(x{k},2)); 
    if usenans, tmp(:) = nan; end
    tmp(beg2(k):end2(k),:) = x{k}(beg1(k):end1(k),:);
    y{k} = tmp;
  end
case 2
  y   = cell(size(x));
  for k = 1:n
    tmp = zeros(size(x{k},1), nsmp(k)); 
    if usenans, tmp(:) = nan; end
    tmp(:, beg2(k):end2(k)) = x{k}(:,beg1(k):end1(k));
    y{k} = tmp;
  end
otherwise
  error('dimensionality of >2 is not supported');
end
