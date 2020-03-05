function [s, cfg] = ft_statfun_diff_itc(cfg, dat, design)

% FT_STATFUN_DIFF_ITC computes the difference in the inter-trial coherence between
% two conditions. The input data for this test should consist of complex-values
% spectral estimates, e.g. computed using FT_FREQANALYSIS with cfg.method='mtmfft',
% 'wavelet' or 'mtmconvcol'.
%
% The ITC is a measure of phase consistency over trials. By randomlly shuffling the
% trials  between the two consitions and repeatedly computing the ITC difference, you
% can test the significance of the two conditions having a different ITC.
%
% A difference in the number of trials poer condition will affect the ITC, however
% since the number of trials remains the same for each random permutation, this bias
% is reflected in the randomization distribution.
%
% Use this function by calling the high-level statistic functions as
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
% with the following configuration option:
%   cfg.statistic = 'ft_statfun_diff_itc'
%
% For this specific statistic there is no known parametric distribution, hence the
% probability and critical value cannot be computed analytically. This specific
% statistic can therefore only be used with cfg.method='montecarlo'. If you want to
% do this in combination with cfg.correctm='cluster', you also need to specify
% cfg.clusterthreshold='nonparametric_common' or 'nonparametric_individual'.
%
% You can specify the following configuration options:
%   cfg.complex = string, 'diffabs' (default) to compute the difference of the absolute ITC values,
%                 or 'absdiff' to compute the absolute value of the difference in the complex ITC values.
%
% The experimental design is specified as:
%   cfg.ivar  = independent variable, row number of the design that contains the labels of the conditions to be compared (default=1)
%
% The labels for the independent variable should be specified as the number 1 and 2.
%
% See also FT_FREQSTATISTICS and FT_STATISTICS_MONTECARLO

% Copyright (C) 2008-2014, Robert Oostenveld
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
cfg.complex        = ft_getopt(cfg, 'complex', 'diffabs');
cfg.ivar           = ft_getopt(cfg, 'ivar', 1);

sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);
df1  = length(sel1);
df2  = length(sel2);

if (df1+df2)<size(design, 2)
  ft_error('inappropriate design, should contain 1''s and 2''s');
end

if isreal(dat)
  ft_error('the data should be complex, i.e. computed with FT_FREQANALYSIS and cfg.output="fourier"');
end

% normalise the complex data in each trial
dat = dat./abs(dat);

switch cfg.complex
  case 'diffabs'
    % first compute the absolute, then take the difference
    % this is not sensitive to phase differences
    itc1 = abs(mean(dat(:,sel1), 2)); % ITC is the length of the average complex number
    itc2 = abs(mean(dat(:,sel2), 2)); % ITC is the length of the average complex number
    s.stat = itc1 - itc2;
  case 'absdiff'
    % first compute the difference, then take the absolute
    % this is sensitive to phase differences
    itc1 = mean(dat(:,sel1), 2); % ITC is here the average complex number
    itc2 = mean(dat(:,sel2), 2); % ITC is here the average complex number
    s.stat = abs(itc1 - itc2);
  otherwise
    ft_error('incorrect specification of cfg.complex');
end
