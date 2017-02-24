function x = cellrowassign(x, y, indx)

for k = 1:numel(x)
  x{k}(indx,:) = y{k};
end 
