function stat = ft_statfun_cohensd(cfg, dat, design)

% FT_STATFUN_COHENSD computes the effect size according to Cohen's d. This function
% supports both unpaired and paired designs.
%
% The table below contains descriptors for magnitudes of Cohen's d.
%   Very small  0.01
%   Small       0.20
%   Medium      0.50
%   Large       0.80
%   Very large  1.20
%   Huge        2.00
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'ft_statfun_cohensd'
%
% The experimental design is specified as:
%   cfg.ivar  = independent variable, row number of the design that contains the labels of the conditions to be compared (default=1)
%   cfg.uvar  = optional, row number of design that contains the labels of the units-of-observation, i.e. subjects or trials (default=2)
%
% The labels for the independent variable should be specified as the number 1 and 2.
% The labels for the unit of observation should be integers ranging from 1 to the
% total number of observations (subjects or trials).
%
% The cfg.uvar option is only needed for paired data, you should leave it empty
% for non-paired data.
%
% See https://en.wikipedia.org/wiki/Effect_size#Cohen.27s_d for a description
% and https://www.psychometrica.de/effect_size.html for an online computation tool.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS

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

% set the defaults
cfg.ivar           = ft_getopt(cfg, 'ivar', 1);
cfg.uvar           = ft_getopt(cfg, 'uvar', []); % default is empty, which means not paired

% from the STD function: a weight of 0 normalizes by N-1 and a weight of 1 normalizes by N
w = 0;

if isempty(cfg.uvar)
  sel1 = find(design(cfg.ivar,:)==1); % select replications that belong to condition 1
  sel2 = find(design(cfg.ivar,:)==2); % select replications that belong to condition 2
  n1 = length(sel1);
  n2 = length(sel2);

  x1 = dat(:,sel1);
  x2 = dat(:,sel2);

  pooled_sd = sqrt(((n1-1)*std(x1, w, 2).^2 + (n2-1)*std(x2, w, 2).^2) / (n1+n2-1));
  cohensd   = (mean(x1, 2) - mean(x2, 2)) ./ pooled_sd;

else
  subj = unique(design(cfg.uvar,:)); % it can also be paired over trials

  n = length(subj);
  sel1 = nan(size(subj));
  sel2 = nan(size(subj));
  for i=1:n
    sel1(i) = find(design(cfg.uvar,:)==subj(i) & design(cfg.ivar,:)==1);
    sel2(i) = find(design(cfg.uvar,:)==subj(i) & design(cfg.ivar,:)==2);
  end
  x1 = dat(:,sel1);
  x2 = dat(:,sel2);

  cohensd = mean(x1 - x2, 2) ./ std(x1 - x2, w, 2);
end

% do not return a real test statistic or a propability
stat.stat = nan(size(dat,1),1);
stat.prob = nan(size(dat,1),1);

% return Cohen's d for all channel-time-frequency points or for all voxels
stat.cohensd = cohensd;

% also return the difference for all channel-time-frequency points or for all voxels
if isempty(cfg.uvar)
  stat.difference = mean(x1,2)-mean(x2,2);
else
  stat.difference = mean(x1-x2,2);
end
