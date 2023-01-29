function test_ft_freqinterpolate

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqinterpolate
% DATATYPE freq

freq.label = {'1';'2'};
freq.freq  = 1:10;
freq.time  = 1:5;
freq.dimord = 'rpt_chan_freq_time';
freq.powspctrm = randn(100,2,10,5);

cfg = [];
cfg.method = 'nan';
cfg.foilim = [4 6];
freqout = ft_freqinterpolate(cfg, freq);