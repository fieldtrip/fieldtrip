function [sel1, sel2] = match_str(a, b, fullout)

% MATCH_STR looks for matching labels in two lists of strings
% and returns the indices into both the 1st and 2nd list of the matches.
% They will be ordered according to the first input argument.
%
% Use as
%   [sel1, sel2] = match_str(strlist1, strlist2)
%
% The strings can be stored as a char matrix or as an vertical array of
% cells, the matching is done for each row.
%
% When including a 1 as the third input argument, the output lists of
% indices will be expanded to the size of the largest input argument.
% Entries that occur only in one of the two inputs will correspond to a 0
% in the output, in this case. This can be convenient in rare cases if the
% size of the input lists is meaningful.

% Copyright (C) 2000-2021, Robert Oostenveld
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

% ensure that both are cell-arrays
if isempty(a)
  a = {};
elseif ~iscell(a)
  a = cellstr(a);
end
if isempty(b)
  b = {};
elseif ~iscell(b)
  b = cellstr(b);
end

% regardless of what optimizations are implemented, the code should remain
% functionally compatible to the original, which is
% sel1 = [];
% sel2 = [];
% for i=1:length(a)
%   for j=1:length(b)
%     if strcmp(a(i),b(j))
%        sel1 = [sel1; i];
%        sel2 = [sel2; j];
%      end
%    end
% end

% ensure that both are column vectors
a = a(:);
b = b(:);
Na = numel(a);
Nb = numel(b);

% this is a very common use pattern that can be dealt with quickly
if isequal(a, b)
  % the returned vectors must be columns
  sel1 = (1:Na)';
  sel2 = (1:Nb)';
  return
end

% According to the original implementation empty numeric elements are
% allowed, but are not returned as match. This is different to empty string
% elements, which are returned as match.
% See also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1808
empty_a = cellfun(@isnumeric, a) & cellfun(@isempty, a);
empty_b = cellfun(@isnumeric, b) & cellfun(@isempty, b);
% the following allows the unique function to operate normally
a(empty_a) = {''};
b(empty_b) = {''};

% replace all unique strings by a unique number and use the fact that
% numeric comparisons are much faster than string comparisons
[dum1, dum2, c] = unique([a; b]);
a = c(1:Na);
b = c((Na+1):end);

% empty numeric elements should never be returned as a match
a(empty_a) = nan;
b(empty_b) = nan;

if nargin < 3 || ~fullout
  sel1 = [];
  sel2 = [];
  for i=1:length(a)
    % s = find(strcmp(a(i), b));  % for string comparison
    s = find(a(i)==b);            % for numeric comparison
    sel2 = [sel2; s];
    s(:) = i;
    sel1 = [sel1; s];
  end
else
  sel1 = zeros(max(Na,Nb),1);
  sel2 = zeros(max(Na,Nb),1);
  for i=1:length(a)
    s = find(a(i)==b);
    sel2(s) = s;
    sel1(s) = i;
  end
end
