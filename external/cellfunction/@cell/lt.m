function [y] = lt(x, val)

y = cell(size(x));
for k = 1:numel(x)
  y{k} = x{k}<val;
end 
