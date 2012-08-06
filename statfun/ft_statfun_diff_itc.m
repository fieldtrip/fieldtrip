function [s] = statfun_diff_itc(cfg, dat, design)

% STATFUN_diff_itc computes the difference in the inter-trial
% coherence (ITC) between two conditions. The input data for this test
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
%   cfg.statistic = 'diff_itc'
%   cfg.method    = 'montecarlo'
% and optionally (for use in the statfun) the option
%  cfg.complex = 'diffabs' to compute the difference of the absolute ITC values, or
%  cfg.complex = 'absdiff' to compute the absolute value of the difference in the complex ITC values
% 
% See FT_FREQSTATISTICS and STATISTICS_MONTECARLO for more details

% Copyright (C) 2008, Robert Oostenveld

% set the defaults
if ~isfield(cfg, 'complex'), cfg.complex = 'diffabs';   end

selA = find(design(cfg.ivar,:)==1);
selB = find(design(cfg.ivar,:)==2);
dfA  = length(selA);
dfB  = length(selB);
if (dfA+dfB)<size(design, 2)
  error('inappropriate design, should contain 1''s and 2''s');
end
if isreal(dat)
  error('the data should be complex, i.e. computed with freqanalysis and cfg.output="fourier"');
end
% normalise the complex data in each trial
dat = dat./abs(dat);

switch cfg.complex
case 'diffabs'
  % first compute the absolute, then take the difference
  % this is not sensitive to phase differences
  itcA = abs(mean(dat(:,selA), 2)); % ITC is the length of the average complex numbers
  itcB = abs(mean(dat(:,selB), 2)); % ITC is the length of the average complex numbers
  s = itcA - itcB;
case 'absdiff'
  % first compute the difference, then take the absolute
  % this is sensitive to phase differences
  itcA = mean(dat(:,selA), 2); % ITC is here the average complex number
  itcB = mean(dat(:,selB), 2); % ITC is here the average complex number
  s = abs(itcA - itcB);
otherwise
  error('incorrect specification of cfg.complex');
end

s.stat = s;

