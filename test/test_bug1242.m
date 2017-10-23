function test_bug1242

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_databrowser

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1242.mat'));

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, timelock);

