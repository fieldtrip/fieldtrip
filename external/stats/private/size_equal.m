function eq = size_equal(varargin)
% returns true if all input arguments are equal to each other
n=numel(varargin);

if n<=1
  eq=true;
  return;
end

size_first=size(varargin{1});

for k=2:n
  if ~isequal(size_first, size(varargin{k}))
    eq=false;
    return;
  end
end

eq=true;
