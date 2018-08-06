function [sel1, sel2] = match_str(a, b)

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
% $Log: match_str.m,v $
% Revision 1.6  2006/11/06 21:11:45  roboos
% also deal with empty [] input
%
% Revision 1.5  2004/11/10 17:11:40  roboos
% reverted to original implementation and reimplemented the speed up
% from scratch. The previous two revisions both were incompatible
% with the original implementation.
%
% Revision 1.4  2004/11/09 15:28:57  roboos
% fixed incompatibility that was introduced by previous speed increase:
% the original version gave back double occurences, and other fieldtrip
% functions (sourceanalysis) rely on this. The previously commited
% version only gave back one occurence of each hit, this is fixed by jansch
% in this version
%
% Revision 1.3  2004/10/22 15:59:41  roboos
% large speed increase by replacing 2 nested for loops by a standard matlab function (intersect)
%
% Revision 1.2  2003/03/17 10:37:28  roberto
% improved general help comments and added copyrights

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

% ensure that both are column vectors
a = a(:);
b = b(:);

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

% replace all unique strings by a unique number and use the fact that
% numeric comparisons are much faster than string comparisons
[dum1, dum2, c] = unique([a; b]);
a = c(1:length(a));
b = c((length(a)+1):end);

sel1 = [];
sel2 = [];
for i=1:length(a)
  % s = find(strcmp(a(i), b));  % for string comparison
  s = find(a(i)==b);            % for numeric comparison
  sel1 = [sel1; repmat(i, size(s))];
  sel2 = [sel2; s];
end
