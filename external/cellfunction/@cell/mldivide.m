function z = mldivide(x,y)

if iscell(x) && iscell(y)
  z = cell(size(x));
  for k = 1:numel(x)
    z{k} = x{k}\y{k};
  end
elseif iscell(x)
  z = cell(size(x));
  for k = 1:numel(x)
    z{k} = x{k}\y;
  end
elseif iscell(y)
  z = cell(size(y));
  for k = 1:numel(y)
    z{k} = x\y{k};
  end
end
