function test_ft_freqdescriptives

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqdescriptives
% DATATYPE freq

freq.label = {'1';'2'};
freq.freq  = 1:10;
freq.time  = 1:5;
freq.dimord = 'rpt_chan_freq_time';
freq.powspctrm = randn(100,2,10,5);

cfg = [];
cfg.variance = 'yes';
freqout = ft_freqdescriptives(cfg, freq);