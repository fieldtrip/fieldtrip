function writelines(lines, filename)

assert(nargin==2); 
assert(ischar(filename));
assert(all(cellfun(@ischar, lines)));

fid = fopen(filename, 'w');
if fid<0
  error('cannot open file');
end

for i=1:length(lines)
  fwrite(fid, lines{i}, 'char');
  fwrite(fid, newline, 'char'); % Go to the next line
end

fclose(fid);