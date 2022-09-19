function test_bug1242

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_databrowser

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1242.mat'));

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, timelock);

