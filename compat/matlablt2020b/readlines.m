function lines = readlines(filename)

% drop in replacement of the readlines function, without additional
% options, i.e. it uses all defaults as compared to the readlines function

newline = {'\r' '\n' '\r\n'};

% Open file in binary mode
fid = fopen(filename, 'rb');
if fid == -1
  error('Could not open file: %s', filename);
end

% Read first few bytes to detect Byte-Order-Marker
firstBytes = fread(fid, 4, '*uint8')';
fseek(fid, 0, 'bof'); % Rewind to start

% Determine encoding based on BOM
if length(firstBytes) >= 3 && isequal(firstBytes(1:3), [239 187 191])
  % UTF-8 BOM
  encoding = 'UTF-8';
  bomLength = 3;
elseif length(firstBytes) >= 2 && isequal(firstBytes(1:2), [255 254])
  % UTF-16LE BOM
  encoding = 'UTF-16LE';
  bomLength = 2;
elseif length(firstBytes) >= 2 && isequal(firstBytes(1:2), [254 255])
  % UTF-16BE BOM
  encoding = 'UTF-16BE';
  bomLength = 2;
elseif length(firstBytes) >= 4 && isequal(firstBytes(1:4), [255 254 0 0])
  % UTF-32LE BOM
  encoding = 'UTF-32LE';
  bomLength = 4;
elseif length(firstBytes) >= 4 && isequal(firstBytes(1:4), [0 0 254 255])
  % UTF-32BE BOM
  encoding = 'UTF-32BE';
  bomLength = 4;
else
  % No BOM, assume UTF-8 (or system default)
  encoding = 'UTF-8';
  bomLength = 0;
end

% Skip BOM if present
if bomLength > 0
  fseek(fid, bomLength, 'bof');
end

% Read the rest of the file as text
text = fread(fid, Inf, '*char')';
fclose(fid);

% Convert to correct encoding
if strcmp(encoding, 'UTF-8')
  lines = strsplit(text, newline);
else
  % For UTF-16/32, use native2unicode if needed
  lines = strsplit(text, newline);
end
lines = string(lines);
