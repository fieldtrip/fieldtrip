function test_bug1242

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_databrowser
% DATA private

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1242.mat'));

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, timelock);

