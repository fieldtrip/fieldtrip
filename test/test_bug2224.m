function test_bug2224

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_selectdata ft_selectdata_new ft_postamble ft_postamble_previous ft_postamble_history

freq1 = [];
freq1.dimord = 'chan_freq';
freq1.freq = 1:10;
freq1.label = {'1', '2', '3'};
freq1.powspctrm = randn(3,10);
freq1.powspctrm(1) = 1;
freq1.cfg.v = 1;

freq2 = [];
freq2.dimord = 'chan_freq';
freq2.freq = 1:10;
freq2.label = {'1', '2', '3'};
freq2.powspctrm = randn(3,10);
freq2.powspctrm(1) = 2;
freq2.cfg.v = 2;

[f1, f2] = ft_selectdata([], freq1, freq2);

assert(isfield(f1.cfg, 'previous'), 'one failed');
assert(isfield(f2.cfg, 'previous'), 'two failed');
assert(numel(f1.cfg.previous)==1,   'three failed');
assert(numel(f2.cfg.previous)==1,   'four failed');
