function [y] = cellsubselect(x, boolvec, dim)

% [Y] = CELLSUBSELECT(X, BOOLVEC, DIM) outputs a cell-arry Y with the same dimension as X
% but from each input cell a subset of rows or columns are kept according
% to boolvec (which assumes the input data to be concatenated across cells)
% 
% X should be a linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells and the sum of samples in
% the other dimension should equal the length of boolvec

scx1 = cellfun('size', x, 1);
scx2 = cellfun('size', x, 2);

if nargin<3
  % reconstruct dim from the boolvec
  if length(boolvec)==sum(scx1), 
    dim  = 1;
  elseif length(boolvec)==sum(scx2),
    dim  = 2;
  elseif all(scx1==scx1(1)) && numel(boolvec==scx1(1))
    dim  = 1;
  elseif all(scx2==scx2(1)) && numel(boolvec==scx2(1))
    dim  = 2;
  else
    error('the length of boolvec should correspond to the summed number of samples across on of the dimensions of the input cells');
  end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellsubselect');
end

switch dim
  case 1
    nsmp = scx1;
  case 2
    nsmp = scx2;
  otherwise
    error('dim argument should be 1 or 2');
end

if numel(boolvec)==sum(nsmp)
  sel = mat2cell(boolvec,1,nsmp(:)');
elseif numel(boolvec)==nsmp(1)
  sel = repmat({boolvec}, nx);
end

dim = repmat({dim}, nx);
y   = cellfun(@subc, x, sel, dim, 'UniformOutput', 0);

function [y] = subc(x, boolvec, dim)

switch dim
  case 1
    y = x(boolvec, :);
  case 2
    y = x(:, boolvec);
otherwise
end
