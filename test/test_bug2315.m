function test_bug2315

% WALLTIME 00:03:00

% TEST test_bug2315
% TEST ft_databrowser ft_prepare_layout


load dataFIC

cfg = [];
cfg.viewmode = '3column';
cfg.channel = 1:15;
ft_databrowser(cfg, dataFIC)

cfg = [];
cfg.viewmode = 'CTF151';
cfg.channel = 'MEG';
ft_databrowser(cfg, dataFIC);

