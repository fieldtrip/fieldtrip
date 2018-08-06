function test_bug2051

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_math

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2051/source_coh_lft.mat'))

cfg = [];
cfg.parameter = 'avg.pow';
cfg.operation = 'log10';
powlog = ft_math(cfg, source_coh_lft);

% sanity check on some other data
timelock = [];
timelock.pow = randn(1,100).^2;
timelock.label = {'a'};
timelock.time = 1:100;
timelock.dimord = 'chan_time';

cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'log10';
powlog = ft_math(cfg, timelock);
