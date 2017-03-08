function test_bug931

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_appendfreq

freq1.label = {'1'};
freq1.time = [1 2];
freq1.freq = [1 2 3];
freq1.dimord = 'chan_freq_time';
freq1.powspctrm = randn(1,3,2);

freq2.label = {'1'};
freq2.time = [1 2];
freq2.freq = [1 2 3] + 1; % shifted by one
freq2.dimord = 'chan_freq_time';
freq2.powspctrm = randn(1,3,2);

cfg = [];
cfg.parameter = 'powspctrm';

try
  % since the time axis are different, it will try to append over the time dimension (and fail)
  ft_appendfreq(cfg, freq1, freq2);
catch
  disp('it produced the expected error');
end



