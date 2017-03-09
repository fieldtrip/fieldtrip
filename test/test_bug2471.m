function test_bug2471

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_timelockgrandaverage

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2471.mat'));

cfg = [];
cfg.channel = 'all';
cfg.keepindividual = 'no';
cfg.method = 'within';
out_data = ft_timelockgrandaverage(cfg,block{:});
