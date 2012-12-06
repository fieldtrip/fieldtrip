function y = subsasgn(x, index, val)

% SUBSASGN Assign a new value to a specified field in a config objects and increment its assignment counter.

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if length(index)==1
  switch index.type
    case '.'
      y = set(x, index.subs, val);
    case '{}'
      error('Cell contents reference from a non-cell array object.');
    case '()'
      error('Index exceeds matrix dimensions.');
    otherwise
      error('Incorrect contents reference');
  end
elseif length(index)>1 && strcmp(index(1).type, '()')
  if ~isfield(x, index(2).subs)
    % an empty field should be added to all elements of the array
    for i=1:numel(x)
      y(i) = subsasgn(x(i), index(2:end), []);
    end
  else
    % the field is already present in the array
    y = x;
  end
  % the value of the field should only be changed for the specific element of the array
  y(index(1).subs{1}) = subsasgn(y(index(1).subs{1}), index(2:end), val);
else
  % use recursion to find the subfield that is being indexed
  if ~isfield(x, index(1).subs)
    % the subfield does not exists, create a new config subfield
    switch index(2).type
      case '.'
        x  = subsasgn(x, index(1), config());
      case '()'
        x  = subsasgn(x, index(1), []);
      case '{}'
        x  = subsasgn(x, index(1), {});
      otherwise
        error('Incorrect contents reference');
    end
  end
  y1 = subsref(x, index(1), false);
  y2 = subsasgn(y1, index(2:end), val);
  y  = set(x, index(1).subs, y2);
end

