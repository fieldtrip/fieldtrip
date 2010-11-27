function [sel1, sel2] = match_str(a, b);

% MATCH_STR looks for matching labels in two listst of strings
% and returns the indices into both the 1st and 2nd list of the matches.
% They will be ordered according to the first input argument.
%
% [sel1, sel2] = match_str(strlist1, strlist2)
%
% The strings can be stored as a char matrix or as an vertical array of
% cells, the matching is done for each row.

% Copyright (C) 2000, Robert Oostenveld
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
% $Id: match_str.m 951 2010-04-21 18:24:01Z roboos $

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
% for i=1:length(a)
%   for j=1:length(b)
%     if strcmp(a(i),b(j))
%       sel1 = [sel1; i];
%       sel2 = [sel2; j];
%     end
%   end
% end

% ensure that both are column vectors
a = a(:);
b = b(:);
Na = numel(a);
Nb = numel(b);

% replace all unique strings by a unique number and use the fact that
% numeric comparisons are much faster than string comparisons
[dum1, dum2, c] = unique([a; b]);
a = c(1:Na);
b = c((Na+1):end);

sel1 = [];
sel2 = [];
for i=1:length(a)
  % s = find(strcmp(a(i), b));  % for string comparison
  s = find(a(i)==b);            % for numeric comparison
  sel2 = [sel2; s];
  s(:) = i;
  sel1 = [sel1; s];
end
