function write_off(filename, pnt, plc)

% WRITE_OFF writes a set of geometrical planar forms (called piecewise linear complex, PLC)
% to an ascii *.off file, which is a file format created by Princeton Shape Benchmark
%
% Use as
%   write_stl(filename, pnt, tri)
%
% See also READ_OFF

% Copyright (C) 2010, Cristiano Micheli
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

nedges = 0;
fid  = fopen(filename, 'wb');
npnt = size(pnt,1);
nplc = size(plc,1);

% check that the indexes of plc are correct (0 convention)
if ~sum(any(plc==0))
  plc = plc - 1;
end

fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d %d\n',npnt,nplc,nedges);
for i=1:npnt
  fprintf(fid, '%f %f %f\n', pnt(i,1), pnt(i,2), pnt(i,3)); 
end
for i=1:nplc
  str = '%d '; 
  nvert = size(plc(i,:),2);
  str2 = 'nvert,';
  for j=1:nvert
    str  = [str '%d '];
    str2 = [str2 'plc(i,' num2str(j) '),'];
  end
  str = [str '\n''']; str2 = str2(1:end-1);
  eval(['fprintf(fid,''' str ',' str2 ');']); 
end

fclose(fid);
