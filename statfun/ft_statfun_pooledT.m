function [s, cfg] = ft_statfun_pooledT(cfg, dat, design)

% FT_STATFUN_POOLEDT computes the pooled t-value over a number of replications. The
% idea is that you compute a contrast between two conditions per subject The t-values
% are pooled over subjects and compared against the pooled pseudo-values. Since
% according to H0 the expected t-value for each subject value is zero, the difference
% between the pooled t-value and the pseudo-value (which is set to zero) is a
% fixed-effects statistic.
% 
% The computation of the difference between pooled t-values can be repeated after
% randomly permuting the t-values and pseudo-values within the subjects. Each random
% permutation gives you an estimate of the difference. The random permutations build
% up a randomization distributin, against which you can compare the observed pooled
% t-values.
% 
% The statistical inference based on the comparison of the observed pooled t-values
% with the randomization distribution is not a fixed-effect statistic, one or a few
% outlier will cause the randomization distribution to broaden and result in the
% conclusion of "not significant".
% 
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option
%   cfg.statistic = 'ft_statfun_pooledT'
%
% Configuration options that are relevant for this function are
%   cfg.ivar      = number, index into the design matrix with the independent variable
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS

% Copyright (C) 2007, Robert Oostenveld
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


selA = find(design(cfg.ivar,:)==1);
selB = find(design(cfg.ivar,:)==2);
dfA  = length(selA);
dfB  = length(selB);
if (dfA+dfB)<size(design, 2)
  ft_error('inappropriate design, should only contain 1''s and 2''s');
end
sumA = sum(dat(:,selA), 2);
sumB = sum(dat(:,selB), 2);
s.stat = (sumA - sumB)./sqrt(dfA);

