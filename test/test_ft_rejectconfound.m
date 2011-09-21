function test_ft_rejectconfound

% TEST test_ft_rejectconfound
% TEST ft_rejectconfound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq1 = [];
freq1.label = {'1' '2'};
freq1.freq  = 1:10;
freq1.dimord = 'rpt_chan_freq';
freq1.powspctrm = randn(20,2,10);

cfg = [];
cfg.confound = randn(20,3);
freq1_out = ft_rejectconfound(cfg, freq1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq2 = [];
freq2.label = {'1' '2'};
freq2.freq  = 1:10;
freq2.time  = 1:5;
freq2.dimord = 'rpt_chan_freq_time';
freq2.powspctrm = randn(20,2,10,5);

cfg = [];
cfg.confound = randn(20,3);
freq2_out = ft_rejectconfound(cfg, freq2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock = [];
timelock.label = {'1' '2'};
timelock.time  = 1:5;
timelock.dimord = 'rpt_chan_time';
timelock.trial = randn(20,2,5);

cfg = [];
cfg.confound = randn(20,3);
timelock_out = ft_rejectconfound(cfg, timelock);

