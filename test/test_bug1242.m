function test_bug1242

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_databrowser
% DATA private

load(dccnpath('/project/3031000.02/test/bug1242.mat'));

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, timelock);

