function [s] = statfun_diff(cfg, dat, design)

% STATFUN_diff computes the difference of the mean in two conditions.
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
% See also STATFUN_MEAN for an other example function

selA = find(design(cfg.ivar,:)==1); % selecton condition 1 or A
selB = find(design(cfg.ivar,:)==2); % selecton condition 2 or B
dfA  = length(selA);
dfB  = length(selB);
if (dfA+dfB)<size(design, 2)
  % there are apparently replications that belong neither to condition 1, nor to condition 2
  warning('inappropriate design, it should only contain 1''s and 2''s');
end
% compute the averages and the difference
avgA = mean(dat(:,selA), 2);
avgB = mean(dat(:,selB), 2);
s = avgA - avgB;

% the stat field is used in STATISTICS_MONTECARLO to make the
% randomization distribution, but you can also return other fields
% which will be passed on to the command line in the end.

s.stat = s;

