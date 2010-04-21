function [pnt, tri, nrm] = read_stl(filename);

% READ_STL reads a triangulation from an ascii *.stl file, which is a file
% format native to the stereolithography CAD software created by 3D Systems.
%
% Use as
%   [pnt, tri, nrm] = read_stl(filename)
%
% See also WRITE_STL

% Copyright (C) 2006, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% solid   testsphere
%   facet normal -0.13 -0.13 -0.98
%     outer loop
%       vertex 1.50000 1.50000 0.00000
%       vertex 1.50000 1.11177 0.05111
%       vertex 1.11177 1.50000 0.05111
%     endloop
%   endfacet
%   ...

fid = fopen(filename, 'rt');

ntri = 0;
while ~feof(fid)
  line = fgetl(fid);
  ntri = ntri + ~isempty(findstr('facet normal', line));
end
fseek(fid, 0, 'bof');

tri = zeros(ntri,3);
pnt = zeros(ntri*3,3);

line = fgetl(fid);
name = sscanf(line, 'solid %s');

for i=1:ntri
  line1 = fgetl(fid);
  line2 = fgetl(fid);
  line3 = fgetl(fid);
  line4 = fgetl(fid);
  line5 = fgetl(fid);
  line6 = fgetl(fid);
  line7 = fgetl(fid);
  i1 = (i-1)*3+1;
  i2 = (i-1)*3+2;
  i3 = (i-1)*3+3;
  tri(i,:) = [i1 i2 i3];
  dum = sscanf(strtrim(line1), 'facet normal %f %f %f'); nrm(i,:) = dum(:)';
  dum = sscanf(strtrim(line3), 'vertex %f %f %f'); pnt(i1,:) = dum(:)';
  dum = sscanf(strtrim(line4), 'vertex %f %f %f'); pnt(i2,:) = dum(:)';
  dum = sscanf(strtrim(line5), 'vertex %f %f %f'); pnt(i3,:) = dum(:)';
end

fclose(fid);

