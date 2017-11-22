function [s, cfg] = ft_statfun_diff_itc(cfg, dat, design)

% FT_STATFUN_DIFF_ITC computes the difference in the inter-trial
% coherence between two conditions. The input data for this test
% should consist of complex-values spectral estimates, e.g. computed
% using FT_FREQANALYSIS with method=mtmfft, wavelet or mtmconvcol.
%
% The ITC is a measure of phase consistency over trials. By randomlly
% shuffling the trials  between the two consitions and repeatedly
% computing the ITC difference, you can test the significance of the
% two conditions having a different ITC.
%
% A difference in the number of trials poer condition will affect the
% ITC, however since the number of trials remains the same for each
% random permutation, this bias is reflected in the randomization
% distribution.
%
% Use this function by calling 
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
% with the following configuration options
%   cfg.method    = 'montecarlo'
%   cfg.statistic = 'diff_itc'
% and optionally the options
%  cfg.complex    = 'diffabs' to compute the difference of the absolute ITC values (default), or
%                   'absdiff' to compute the absolute value of the difference in the complex ITC values.
% 
% NOTE: For this specific statistic there is no known parametric distribution, hence
% the probability and critical value cannot be computed. This specific statistic can
% therefore only be used with cfg.method='montecarlo'. If you want to do this in combination
% with cfg.correctm='cluster', you need cfg.clusterthreshold='nonparametric_common' or
% cfg.clusterthreshold='nonparametric_individual'.
%
% See FT_FREQSTATISTICS and FT_STATISTICS_MONTECARLO for more details

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
if ~isfield(cfg, 'complex'), cfg.complex = 'diffabs';   end

selA = find(design(cfg.ivar,:)==1);
selB = find(design(cfg.ivar,:)==2);
dfA  = length(selA);
dfB  = length(selB);
if (dfA+dfB)<size(design, 2)
  ft_error('inappropriate design, should contain 1''s and 2''s');
end
if isreal(dat)
  ft_error('the data should be complex, i.e. computed with freqanalysis and cfg.output="fourier"');
end
% normalise the complex data in each trial
dat = dat./abs(dat);

switch cfg.complex
case 'diffabs'
  % first compute the absolute, then take the difference
  % this is not sensitive to phase differences
  itcA = abs(mean(dat(:,selA), 2)); % ITC is the length of the average complex numbers
  itcB = abs(mean(dat(:,selB), 2)); % ITC is the length of the average complex numbers
  s.stat = itcA - itcB;
case 'absdiff'
  % first compute the difference, then take the absolute
  % this is sensitive to phase differences
  itcA = mean(dat(:,selA), 2); % ITC is here the average complex number
  itcB = mean(dat(:,selB), 2); % ITC is here the average complex number
  s.stat = abs(itcA - itcB);
otherwise
  ft_error('incorrect specification of cfg.complex');
end

