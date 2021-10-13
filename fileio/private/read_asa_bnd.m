function bnd = read_asa_bnd(fn)

% READ_ASA_BND reads an ASA boundary triangulation file
% converting the units of the vertices to mm

% Copyright (C) 2002, Robert Oostenveld
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

Npnt = read_ini(fn, 'NumberPositions=', '%d');
Ntri = read_ini(fn, 'NumberPolygons=', '%d');
Unit = read_ini(fn, 'UnitPosition', '%s');

pnt = read_ini(fn, 'Positions', '%f');
if any(size(pnt)~=[Npnt,3])
  pnt_file = read_ini(fn, 'Positions', '%s');
  [path, name, ext] = fileparts(fn);
  fid = fopen_or_error(fullfile(path, pnt_file), 'rb', 'ieee-le');
  pnt = fread(fid, [3,Npnt], 'float')';
  fclose(fid);
end

tri = read_ini(fn, 'Polygons', '%f');
if any(size(tri)~=[Ntri,3])
  tri_file = read_ini(fn, 'Polygons', '%s');
  [path, name, ext] = fileparts(fn);
  fid = fopen_or_error(fullfile(path, tri_file), 'rb', 'ieee-le');
  tri = fread(fid, [3,Ntri], 'int32')';
  fclose(fid);
end

if strcmpi(Unit,'mm')
  pnt   = 1*pnt;
elseif strcmpi(Unit,'cm')
  pnt   = 100*pnt;
elseif strcmpi(Unit,'m')
  pnt   = 1000*pnt;
else
  ft_error(sprintf('Unknown unit of distance for triangulated boundary (%s)', Unit));
end

bnd.pnt = pnt;
bnd.tri = tri + 1;      % 1-offset instead of 0-offset
% bnd.tri = fliplr(bnd.tri);    % invert order of triangle vertices

