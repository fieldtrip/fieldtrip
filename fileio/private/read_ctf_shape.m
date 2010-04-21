function [shape] = read_ctf_shape(filename);

% READ_CTF_SHAPE reads headshape points and header information
% from a CTF *.shape teh accompanying *.shape_info file.
%
% Use as
%   [shape] = read_ctf_shape(filename)
% where filename should have the .shape extension 

% Copyright (C) 2003, Robert Oostenveld
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

shape = read_ctf_ascii([filename '_info']);

if ~strcmp(shape.MRI_Info.COORDINATES, 'HEAD')
  warning('points on head shape are NOT in headcoordinates')
end

fid = fopen(filename, 'rt');
num = fscanf(fid, '%d', 1);
shape.pnt = fscanf(fid, '%f', inf);
shape.pnt = reshape(shape.pnt, [3 num])';
fclose(fid);

