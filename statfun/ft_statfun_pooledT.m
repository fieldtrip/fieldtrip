function [s, cfg] = ft_statfun_pooledT(cfg, dat, design)

% FT_STATFUN_POOLEDT computes the pooled t-value over a number of replications. The
% idea behind this function is that you first (prior to calling this function)
% compute a contrast between two conditions per subject, and that subsequently you
% test this over subjects using random sign-flipping.
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option
%   cfg.statistic = 'ft_statfun_pooledT'
%
% The expected values for the pooled-t, which is zero according to H0, have to be
% passed as pseudo-values. The subject-specific t-values will be randomly swapped with
% the pseudo-values and the difference is computed; in effect this implements random
% sign-flipping.
%
% The randimization distribution (with optional clustering) of the randomly
% sign-flipped pooled-t values is computed and used for statistical inference.
%
% Note that, although the output of this function is to be interpreted as a
% fixed-effects statistic, the statistical inference based on the comparison of the
% observed pooled t-values with the randomization distribution is not a fixed-effect
% statistic, one or a few outlier will cause the randomization distribution to
% broaden and result in the conclusion of "not significant".
%
% The experimental design is specified as:
%   cfg.ivar  = independent variable, row number of the design that contains the labels of the conditions to be sign-flipped (default=1)
%
% The labels independent variable should be specified as the number 1 for the
% observed t-values and 2 for the pseudo-values.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS

% Copyright (C) 2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org for the
% documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify it under the
%    terms of the GNU General Public License as published by the Free Software
%    Foundation, either version 3 of the License, or (at your option) any later
%    version.
%
%    FieldTrip is distributed in the hope that it will be useful, but WITHOUT ANY
%    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%    PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License along with
%    FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% set the defaults
cfg.ivar = ft_getopt(cfg, 'ivar', 1);

sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);
df1  = length(sel1);
df2  = length(sel2);

if (df1+df2)<size(design, 2)
  ft_error('inappropriate design, it should only contain 1''s and 2''s');
end

if df1~=df2
  ft_error('inappropriate design, it should contain the same number of 1''s and 2''s');
end

sum1 = sum(dat(:,sel1), 2);
sum2 = sum(dat(:,sel2), 2);

% return the pooled-t value as the statistic of interest
s.stat = (sum1 - sum2)./sqrt(df1);
