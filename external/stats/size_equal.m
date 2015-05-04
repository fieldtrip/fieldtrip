function eq = size_equal(varargin)
n=numel(varargin);

if n<=1
    eq=true;
    return;
end

size_first=size(varargin{1});

for k=2:n
    if ~isequal(size_first, varargin{k})
        eq=false;
        return;
    end
end

eq=true;

end
