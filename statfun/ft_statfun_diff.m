function [s, cfg] = ft_statfun_diff(cfg, dat, design)

% FT_STATFUN_DIFF demonstrates how to compute the difference of the mean in two
% conditions. Although it can be used for statistical testing, it will have rather
% limited sensitivity and is not really suited for inferential testing.
%
% This function serves as an example for a statfun. You can use such a function with
% the statistical framework in FieldTrip using FT_TIMELOCKSTATISTICS,
% FT_FREQSTATISTICS or FT_SOURCESTATISTICS to perform a statistical test, without
% having to worry about the representation of the data.
%
% Use this function by calling the high-level statistic functions as
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
% with the following configuration option:
%   cfg.statistic = 'ft_statfun_diff_itc'
%
% The experimental design is specified as:
%   cfg.ivar  = independent variable, row number of the design that contains the labels of the conditions to be compared (default=1)
%
% The labels for the independent variable should be specified as the number 1 and 2.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS, and see FT_STATFUN_MEAN for a similar example

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

% set the defaults
cfg.ivar           = ft_getopt(cfg, 'ivar', 1);

sel1 = find(design(cfg.ivar,:)==1); % select replications that belong to condition 1
sel2 = find(design(cfg.ivar,:)==2); % select replications that belong to condition 2
n1  = length(sel1);
n2  = length(sel2);

if (n1+n2)<size(design, 2)
  % there are apparently replications that belong neither to condition 1, nor to condition 2
  ft_warning('inappropriate design, it should only contain 1''s and 2''s');
end

% compute the averages and the difference
avg1 = nanmean(dat(:,sel1), 2);
avg2 = nanmean(dat(:,sel2), 2);

% return the difference in the means as statistic of interest
s.stat = avg1 - avg2;
