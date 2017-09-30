function [s, cfg] = ft_statfun_diff(cfg, dat, design)

% FT_STATFUN_DIFF computes the difference of the mean in two conditions.
% Although it can be used for statistical testing, it is not very
% usefull since it will have rather limited sensitivity.
% 
% The purpose of this function is to show you an example on how to
% write a statfun that expresses the difference in the data between
% two conditions. You can use such a function with the statistical
% framework in FieldTrip to perform a simple (or more complex)
% permutation test, without having to worry about the representation
% of the data.
%
% See also FT_STATFUN_MEAN for another example function

% Copyright (C) 2006, Robert Oostenveld 
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

selA = find(design(cfg.ivar,:)==1); % selecton condition 1 or A
selB = find(design(cfg.ivar,:)==2); % selecton condition 2 or B
dfA  = length(selA);
dfB  = length(selB);
if (dfA+dfB)<size(design, 2)
  % there are apparently replications that belong neither to condition 1, nor to condition 2
  ft_warning('inappropriate design, it should only contain 1''s and 2''s');
end
% compute the averages and the difference
avgA = nanmean(dat(:,selA), 2);
avgB = nanmean(dat(:,selB), 2);
s.stat = avgA - avgB;

% the stat field is used in STATISTICS_MONTECARLO to make the
% randomization distribution, but you can also return other fields
% which will be passed on to the command line in the end.

