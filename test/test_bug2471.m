function test_bug2471

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_timelockgrandaverage
% DATA private

load(dccnpath('/project/3031000.02/test/bug2471.mat'));

cfg = [];
cfg.channel = 'all';
cfg.keepindividual = 'no';
cfg.method = 'within';
out_data = ft_timelockgrandaverage(cfg,block{:});
