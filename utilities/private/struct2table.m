function table = struct2table(s)

% STRUCT2TABLE converts a struct-array to a cell-array of strings that represents a table
%
% Example
%   s(1).a = 1
%   s(1).b = 2
%   s(2).a = 3
%   s(2).b = 4
%   disp(struct2table(s))
%
% See also PRINTSTRUCT, APPENDSTRUCT

% Copyright (C) 2017, Robert Oostenveld
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

fn = fieldnames(s);
colwidth = zeros(size(fn));

% first ensure that all fields are strings
for i=1:numel(s)
  for j=1:numel(fn)
    val = s(i).(fn{j});
    switch class(val)
      case 'char'
        % nothing to be done
      case {'single', 'double'}
        val = num2str(val);
      case 'logical'
        if val
          val = 'true';
        else
          val = 'false';
        end
      otherwise
        error('not yet implemented for class "%s"', class(val));
    end
    s(i).(fn{j}) = val;
  end
end

for i=1:numel(fn)
  colwidth(i) = max(cellfun(@length, {s.(fn{i})}));
end

header = '|';
for i=1:numel(fn)
  colname = fn{i};
  if numel(colname)>1 && colname(1)=='X' && colname(end)=='X'
    % the name of the field has been base64 encoded
    % we want to use the original name for the column header
    colname = fixname(fn{i}, 'X_base64decode_X');
  end
  % update the width of the column to the header
  colwidth(i) = max(colwidth(i), length(colname));
  header = [header pad(colname, colwidth(i)+1, 'left', ' '), ' |'];
end

divider = repmat('-', size(header));

% divider = header;
% divider(divider~='|') = '-';

line = cell(numel(s),1);
for i=1:numel(s)
  line{i} = '|';
  for j=1:numel(fn)
    line{i} = [line{i} pad(s(i).(fn{j}),colwidth(j)+1,'left',' '), ' |'];
  end
end

table = cat(1, header, divider, line);
