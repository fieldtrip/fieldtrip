function [x] = ft_struct2char(x, maxdepth)

% FT_STRUCT2CHAR converts all string elements in a structure
% into char-arrays.
%
% Use as
%   x = ft_struct2char(x)
%
% See also FT_STRUCT2STRING, FT_STRUCT2SINGLE, FT_STRUCT2DOUBLE

% Copyright (C) 2018, Robert Oostenveld
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
  % only convert up to the specified level
  return
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
    
  case 'string'
    % convert the string into a char-array
    a = char(a);
    
  case 'char'
    % keep as it is
    
  otherwise
    % do nothing
end
