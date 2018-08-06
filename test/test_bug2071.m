function test_bug2071

% WALLTIME 00:10:00
% MEM 1gb


% TEST ft_postamble
% TEST ft_postamble_history

data1 = [];
data1.label = {'1', '2', '3', '4'};
data1.powspctrm = rand(4, 10, 10);
data1.freq = 1:10;
data1.time = .1:.1:1;
data1.dimord = 'chan_freq_time';
data1.cfg.this = '1a';
data1.cfg.previous = '1b';

data2 = [];
data2.label = {'2', '3', '1', '4'}; % different order
data2.powspctrm = rand(4, 10, 10);
data2.freq = 1:10;
data2.time = .1:.1:1;
data2.dimord = 'chan_freq_time';
data2.cfg.this = '2a';
data2.cfg.previous = '2b';

data3 = [];
data3.label = {'3', '4', '1', '2'}; % channel missing
data3.powspctrm = rand(4, 10, 10);
data3.freq = 1:10;
data3.time = .1:.1:1;
data3.dimord = 'chan_freq_time';
data3.cfg.this = '3a';
data3.cfg.previous = '3b';

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'add';
dataout = ft_math(cfg, data1, data2, data3);

ft_postamble history dataout
