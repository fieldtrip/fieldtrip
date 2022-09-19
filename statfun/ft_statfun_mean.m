function [s, cfg] = ft_statfun_mean(cfg, dat, design)

% FT_STATFUN_MEAN demonstrates how to compute the mean over all conditions in the
% data. Since this does NOT depend on the experimental design, it cannot be used for
% testing for differences between conditions.
%
% This function serves as an example for a statfun. You can use such a function with
% the statistical framework in FieldTrip using FT_TIMELOCKSTATISTICS,
% FT_FREQSTATISTICS or FT_SOURCESTATISTICS to perform a statistical test, without
% having to worry about the representation of the data.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS, and see FT_STATFUN_DIFF for a similar example

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if numel(unique(design(cfg.ivar,:)))>1
  ft_warning('inappropriate design, it should only contain one condition');
end

% return the means as the statistic of interest
s.stat = mean(dat,2);
