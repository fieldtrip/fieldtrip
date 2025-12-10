function test_bug2559


% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_databrowser ft_checkdata ft_datatype
% DATA private

load(dccnpath('/project/3031000.02/test/bug2559.mat'));

cfg = [];
cfg.layout = 'easycapM10';
ft_databrowser(cfg, timelock);

