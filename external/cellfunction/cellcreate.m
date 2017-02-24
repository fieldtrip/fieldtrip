function x = cellcreate(s1, s2)

x = cell(size(s1));
for k = 1:numel(x)
  x{k} = nan + zeros(s1(k),s2(k));
end

