function test_bug3256

% WALLTIME 00:10:00
% MEM 1gb


% create some ramdom data
cfg = [];
cfg.numtrl = 30;
raw = ft_freqsimulation(cfg);

cfg = [];
timelock = ft_timelockanalysis(cfg, raw);

%%
cfg = [];
cfg.method = 'pca';

raw_comp = ft_componentanalysis(cfg, raw);
assert(ft_datatype(raw_comp, 'raw+comp'));

timelock_comp = ft_componentanalysis(cfg, timelock);
assert(ft_datatype(timelock_comp, 'timelock+comp'));

%%
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
raw_comp_freq = ft_freqanalysis(cfg, raw_comp);
assert(ft_datatype(raw_comp_freq, 'freq+comp'));

timelock_comp_freq = ft_freqanalysis(cfg, raw_comp);
assert(ft_datatype(timelock_comp_freq, 'freq+comp'));

