function test_bug2315

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_databrowser ft_prepare_layout

load(dccnpath('/home/common/matlab/fieldtrip/data/test/dataFIC.mat'));

%%

cfg = [];
cfg.viewmode = 'ordered';
cfg.columns = 5;
cfg.channel = 1:15;
ft_databrowser(cfg, dataFIC)

%%

cfg = [];
cfg.viewmode = 'CTF151';
cfg.channel = 'MEG';
ft_databrowser(cfg, dataFIC);

