function writelines(lines, filename)

% This is a very limited overloaded version of the MathWorks function WRITELINES
% which has been introduced in MATLAB 2022a. This version does not offer full
% backwards compatibility, it is only capable of writing a cell-array with lines
% to a text file.
%
% The directory containing this function shoudl only be added to the path on
% MATLAB versions prior to 2022a.

% Copyright (C) 2023, Konstantinos Tsilimparis

assert(nargin==2, 'only two input arguments are supported'); 
assert(ischar(filename), 'the filename must be a string');
assert(all(cellfun(@ischar, lines)), 'invalid lines in the cell-array');

fid = fopen(filename, 'w');
if fid<0
  error('cannot open file "%s"', filename);
end

for i=1:length(lines)
  fwrite(fid, lines{i}, 'char');
  fwrite(fid, newline, 'char'); % Go to the next line
end

fclose(fid);
