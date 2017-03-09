function test_bug2509

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_selectdata

freq = [];
freq.dimord = 'rpt_chan_freq_time';
freq.label = {'1', '2', '3'};
freq.freq = 1:4;
freq.time = 1:5;
freq.powspctrm = randn(2, 3, 4, 5);
 
cfg = [];
cfg.avgovertime = 'yes';
output = ft_selectdata(cfg, freq);
assert(isequal(output.time, 3))
