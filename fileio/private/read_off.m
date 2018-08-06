function [pnt, tri] = read_off(fn)

% READ_OFF reads vertices and triangles from a OFF format triangulation file
%
% [pnt, tri] = read_off(filename)
%
% See also READ_TRI, READ_BND

% Copyright (C) 1998, Robert Oostenveld
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

fid = fopen(fn, 'rt');
if fid~=-1

  % scan the file type
  [s, count] = fscanf(fid, '%s\n', 1);
  if ~strcmp(s,'OFF')
    msg = sprintf('wrong file type %s', s);
    ft_error(msg);
  end

  % read the number of vertex points and triangles
  [val, count] = fscanf(fid, '%d', 3);
  Npnt = val(1)
  Ntri = val(2)

  % read the vertex points
  pnt  = fscanf(fid, '%f', [3, Npnt]);
  pnt  = pnt(1:3,:)';

  % read the triangles
  tri = fscanf(fid, '%d', [4, Ntri]);
  tri = (tri(2:4,:)+1)';
  fclose(fid);

else
  ft_error('unable to open file');
end


