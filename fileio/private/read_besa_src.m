function [src] = read_besa_src(filename);

% READ_BESA_SRC reads a beamformer source reconstruction from a BESA file
%
% Use as
%   [src] = read_besa_src(filename)
%
% The output structure contains a minimal representation of the contents
% of the file.

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: read_besa_src.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.2  2005/10/05 06:32:08  roboos
% removed forced reading of (incorrect) condition label
%
% Revision 1.1  2005/09/01 09:22:35  roboos
% new implementation, only tested on a single file
%

src = [];
fid = fopen(filename, 'rt');

% read the header information
line = fgetl(fid);
while isempty(strmatch('BESA_SA_IMAGE', line)), line = fgetl(fid); check_feof(fid, filename); end
% while isempty(strmatch('sd'           , line)), line = fgetl(fid); check_feof(fid, filename); end   % this line contains condition information
while isempty(strmatch('Grid dimen'   , line)), line = fgetl(fid); check_feof(fid, filename); end
while isempty(strmatch('X:'           , line)), line = fgetl(fid); check_feof(fid, filename); end; linex = line;
while isempty(strmatch('Y:'           , line)), line = fgetl(fid); check_feof(fid, filename); end; liney = line;
while isempty(strmatch('Z:'           , line)), line = fgetl(fid); check_feof(fid, filename); end; linez = line;

tok = tokenize(linex, ' '); 
src.X = [str2num(tok{2}) str2num(tok{3}) str2num(tok{4})];
tok = tokenize(liney, ' '); 
src.Y = [str2num(tok{2}) str2num(tok{3}) str2num(tok{4})];
tok = tokenize(linez, ' '); 
src.Z = [str2num(tok{2}) str2num(tok{3}) str2num(tok{4})];

nx = src.X(3);
ny = src.Y(3);
nz = src.Z(3);

src.vol = zeros(nx, ny, nz);

for i=1:nz
  % search up to the next slice
  while isempty(strmatch(sprintf('Z: %d', i-1), line)), line = fgetl(fid); check_feof(fid, filename); end;
  % read all the values for this slice
  buf = fscanf(fid, '%f', [nx ny]);
  src.vol(:,:,i) = buf;
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_feof(fid, filename)
if feof(fid)
  error(sprintf('could not read all information from file ''%s''', filename));
end

