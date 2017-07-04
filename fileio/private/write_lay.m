function writelay(filename, layout)

fid = fopen(filename, 'wt');

for i=1:numel(layout.label)
  fprintf(fid, '%3d %12f %12f %12f %12f %15s\n', i, layout.pos(i,1), layout.pos(i,2), layout.width(i), layout.height(i), layout.label{i});
end

fclose(fid)

