function test_bug2315

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_databrowser ft_prepare_layout

load(dccnpath('/home/common/matlab/fieldtrip/data/test/dataFIC.mat'));

cfg = [];
cfg.viewmode = '3column';
cfg.channel = 1:15;
ft_databrowser(cfg, dataFIC)

cfg = [];
cfg.viewmode = 'CTF151';
cfg.channel = 'MEG';
ft_databrowser(cfg, dataFIC);

