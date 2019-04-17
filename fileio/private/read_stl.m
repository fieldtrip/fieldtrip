function [pnt, tri, nrm] = read_stl(filename)

% READ_STL reads a triangulation from an ascii or binary *.stl file, which
% is a file format native to the stereolithography CAD software created by
% 3D Systems.
%
% Use as
%   [pnt, tri, nrm] = read_stl(filename)
%
% The format is described at http://en.wikipedia.org/wiki/STL_(file_format)
%
% See also WRITE_STL

% Copyright (C) 2006-2011, Robert Oostenveld
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

fid = fopen(filename, 'rt');

% read a small section to determine whether it is ascii or binary
% a binary STL file has an 80 byte asci header, followed by non-printable characters

section = fread(fid, 160, 'uint8');
fseek(fid, 0, 'bof');

if printableascii(section)
  % the first 160 characters are printable ascii, so assume it is an ascii format

  % solid   testsphere
  %   facet normal -0.13 -0.13 -0.98
  %     outer loop
  %       vertex 1.50000 1.50000 0.00000
  %       vertex 1.50000 1.11177 0.05111
  %       vertex 1.11177 1.50000 0.05111
  %     endloop
  %   endfacet
  %   ...

  ntri = 0;
  while ~feof(fid)
    line = fgetl(fid);
    ntri = ntri + contains('facet normal', line);
  end
  fseek(fid, 0, 'bof');

  tri = zeros(ntri,3);
  nrm = zeros(ntri,3);
  pnt = zeros(ntri*3,3);

  line = fgetl(fid);
  name = sscanf(line, 'solid %s');

  for i=1:ntri
    line1 = fgetl(fid);
    line2 = fgetl(fid); % outer loop
    line3 = fgetl(fid);
    line4 = fgetl(fid);
    line5 = fgetl(fid);
    line6 = fgetl(fid); % endloop
    line7 = fgetl(fid); % endfacet
    i1 = (i-1)*3+1;
    i2 = (i-1)*3+2;
    i3 = (i-1)*3+3;
    tri(i,:) = [i1 i2 i3];
    dum = sscanf(strtrim(line1), 'facet normal %f %f %f'); nrm(i,:) = dum(:)';
    dum = sscanf(strtrim(line3), 'vertex %f %f %f'); pnt(i1,:) = dum(:)';
    dum = sscanf(strtrim(line4), 'vertex %f %f %f'); pnt(i2,:) = dum(:)';
    dum = sscanf(strtrim(line5), 'vertex %f %f %f'); pnt(i3,:) = dum(:)';
  end

else
  % reopen the file in binary mode, which does not make a difference on
  % UNIX but it does on windows
  fclose(fid);
  fid = fopen(filename, 'rb');

  fseek(fid, 80, 'bof'); % skip the ascii header
  ntri = fread(fid, 1, 'uint32');
  tri = reshape(1:(ntri*3),[3 ntri])';
  tmp = fread(fid, [12 ntri], '12*float32', 2); % read 12 floats at a time, and skip 2 bytes.
  nrm = tmp(1:3,:)';

  tmp = reshape(tmp(4:end,:),[3 3 ntri]); % position info
  tmp = permute(tmp,[2 3 1]);
  pnt = reshape(tmp, [], 3);

  % the above replaces the below, which is much slower, because it is using
  % a for loop across triangles
%   tri  = zeros(ntri,3);
%   nrm  = zeros(ntri,3);
%   pnt  = zeros(ntri*3,3);
%   attr = zeros(ntri,1);
%   for i=1:ntri
%     i1 = (i-1)*3+1;
%     i2 = (i-1)*3+2;
%     i3 = (i-1)*3+3;
%     tri(i,:) = [i1 i2 i3];
%     nrm(i,:)  = fread(fid, 3, 'float32');
%     pnt(i1,:) = fread(fid, 3, 'float32');
%     pnt(i2,:) = fread(fid, 3, 'float32');
%     pnt(i3,:) = fread(fid, 3, 'float32');
%     attr(i)   = fread(fid, 1, 'uint16'); % Attribute byte count, don't know what it is
%   end % for each triangle
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function retval = printableascii(num)
% Codes 20hex (32dec) to 7Ehex (126dec), known as the printable characters,
% represent letters, digits, punctuation marks, and a few miscellaneous
% symbols. There are 95 printable characters in total.
num = double(num);
num(num==double(sprintf('\n'))) = double(sprintf(' '));
num(num==double(sprintf('\r'))) = double(sprintf(' '));
num(num==double(sprintf('\t'))) = double(sprintf(' '));
retval = all(num>=32 & num<=126);
