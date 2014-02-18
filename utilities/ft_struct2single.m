function [x] = ft_struct2single(x, maxdepth)

% FT_STRUCT2SINGLE converts all double precision numeric data in a structure
% into single precision, which takes up half the amount of memory compared
% to double precision. It will also convert plain matrices and cell-arrays.
%
% Use as
%    x = ft_struct2single(x);
%
% Starting from MATLAB 7.0, you can use single precision data in your
% computations, i.e. you do not have to convert back to double precision.
%
% MATLAB version 6.5 and older only support single precision for storing
% data in memory or on disk, but do not allow computations on single
% precision data. After reading a single precision structure from file, you
% can convert it back with FT_STRUCT2DOUBLE.
%
% See also FT_STRUCT2DOUBLE

% Copyright (C) 2005-2014, Robert Oostenveld
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

if nargin<2
  maxdepth = inf;
end

% convert the data, work recursively through the complete structure
x = convert(x, 0, maxdepth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this subfunction does the actual work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a] = convert(a, depth, maxdepth)

if depth>maxdepth
  error('recursive depth exceeded');
end

switch class(a)
  case 'struct'
    % process all fields of the structure recursively
    fna = fieldnames(a);
    % process all elements of the array
    for j=1:length(a(:))
      % warning, this is a recursive call to traverse nested structures
      for i=1:length(fna)
        fn = fna{i};
        ra = getfield(a(j), fn);
        ra = convert(ra, depth+1, maxdepth);
        a(j) = setfield(a(j), fn, ra);
      end
    end

  case 'cell'
    % process all elements of the cell-array recursively
    % warning, this is a recursive call to traverse nested structures
    for i=1:length(a(:))
      a{i} = convert(a{i}, depth+1, maxdepth);
    end

  case {'double' 'int64' 'uint64' 'int32' 'uint32' 'int16' 'uint16' 'int8' 'uint8'}
    % convert the values to single precision
    if ~issparse(a)
      a = single(a);
    end

  otherwise
    % do nothing
end
