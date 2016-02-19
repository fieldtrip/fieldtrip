function [src] = read_besa_src(filename)

% READ_BESA_SRC reads a beamformer source reconstruction from a BESA file
%
% Use as
%   [src] = read_besa_src(filename)
%
% The output structure contains a minimal representation of the contents
% of the file.

% Copyright (C) 2005, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

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

