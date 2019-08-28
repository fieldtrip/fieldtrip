function t3 = merge_table(t1, t2, key)

% MERGE_TABLE merges two tables where the rows and columns can be partially
% overlapping or different. Values from the 2nd input have precedence in case the
% same row and column is also present in the 1st.
%
% Use as
%   t3 = merge_table(t1, t2, key)
%
% See also TABLE

% Copyright (C) 2019, Robert Oostenveld
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

assert(nargin==3);

c1 = t1.Properties.VariableNames;
c2 = t2.Properties.VariableNames;

assert(contains(key, c1));
assert(contains(key, c2));

[m1, n1] = size(t1);
[m2, n2] = size(t2);

for i=1:n1
  if ischar(t1{1,i})
    t1.(c1{i}) = cellstr(t1{:,i});
  elseif isnumeric(t1{1,i})
    t1.(c1{i}) = num2cell(t1{:,i});
  end
end

for i=1:n2
  if ischar(t2{1,i})
    t2.(c2{i}) = cellstr(t2{:,i});
  elseif isnumeric(t2{1,i})
    t2.(c2{i}) = num2cell(t2{:,i});
  end
end

t3 = table();

% deal with the columns that are present in both inputs
col = intersect(c1, c2);

for i=1:numel(col)
  % copy the existing column from t1 into t3
  t3.(col{i}) = t1.(col{i});
end

for j=1:m2
  row = findkey(t1.(key), t2.(key){j});
  if ~any(row)
    % add an empty row to t3
    row = size(t3,1)+1;
    try
      % on MATLAB R2019a it is required to add a complete empty row first
      t3(row,:) = cell(1, size(t3,2));
    catch
      % on MATLAB R2016b it is not possible to add an empty row, but looping over columns works fine
    end
  end
  for i=1:numel(col)
    % this also adds the row (when required) for R2016b
    t3(row, col{i}) = t2(j, col{i});
  end
end

n3 = size(t3,1);

% deal with the columns that are only present in the 1st input
col = setdiff(c1, c2);

for i=1:numel(col)
  % add an empty column to t3
  t3.(col{i}) = cell(n3, 1);
end

for j=1:m1
  row = findkey(t3.(key), t1.(key){j});
  for i=1:numel(col)
    t3(row, col{i}) = t1(j, col{i});
  end
end

% deal with the columns that are only present in the 2nd input
col = setdiff(c2, c1);

for i=1:numel(col)
  % add an empty column to t3
  t3.(col{i}) = cell(n3, 1);
end

for j=1:m2
  row = findkey(t3.(key), t2.(key){j});
  for i=1:numel(col)
    t3(row, col{i}) = t2(j, col{i});
  end
end

% reorder the columns to match the input
c3 = unique(cat(2, c1, c2), 'stable');
[c3, ia, ib] = intersect(c3, t3.Properties.VariableNames, 'stable');
t3 = t3(:,ib);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this also allows comparisons between different classes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret = findkey(list, elem)
match = @(x) isequal(x, elem);
ret = cellfun(match, list);
