function test_bug1242

% MEM 1500mb
% WALLTIME 00:03:03

% TEST test_bug1242
% TEST ft_databrowser

load(dccnfilename('/home/common/matlab/fieldtrip/data/test/bug1242.mat'));

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, timelock);

