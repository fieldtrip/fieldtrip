function t = join_str(separator,cells)
if isempty(cells)
  t = '';
  return;
end

if ischar(cells)
  t = cells;
  return;
end

t = char(cells{1});

for i=2:length(cells)
  t = [t separator char(cells{i})];
end