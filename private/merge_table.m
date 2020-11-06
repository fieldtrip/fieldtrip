function t3 = merge_table(t1, t2, key)

% MERGE_TABLE merges two tables where the rows and columns can be partially
% overlapping or different. Values from the 2nd input have precedence in case the
% same row and column is also present in the 1st.
%
% Use as
%   t3 = merge_table(t1, t2)
% or
%   t3 = merge_table(t1, t2, key)
%
% See also TABLE, JOIN, INNERJOIN, OUTERJOIN

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


% deal with the easy cases
if isequal(t1, t2)
  t3 = t1;
  return
elseif isempty(t1)
  t3 = t2;
  return
elseif isempty(t2)
  t3 = t1;
  return
end

if nargin>2
  % see https://www.diffen.com/difference/Inner_Join_vs_Outer_Join
  t3 = outerjoin(t1, t2, 'keys', {key}, 'mergekeys', true);
else
  t3 = vertcat(t1, t2);
end
