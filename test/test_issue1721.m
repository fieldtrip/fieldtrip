function test_issue1721

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_multiplotER ft_singleplotER ft_topoplotER ft_multiplotTFR ft_singleplotTFR ft_topoplotTFR ft_topoplotIC topoplot_common

%%

cfg = [];
cfg.layout = 'elec1020.lay';
layout = ft_prepare_layout(cfg);

%%

timelock1 = [];
timelock1.label = layout.label(1:21); % without COMNT and SCALE
timelock1.time = (1:500)/500;
timelock1.avg = randn(21, 500);
timelock1.dimord = 'chan_time';

timelock2 = [];
timelock2.label = layout.label(1:21); % without COMNT and SCALE
timelock2.time = (1:500)/500;
timelock2.avg = randn(21, 500);

%%

cfg = [];
cfg.layout = layout;
ft_multiplotER(cfg, timelock1, timelock2);

% the bug was:
%  - Select some channels in the multiplot and click. This opens figure 2 and 3, where 3 contains the singleplot.
%  - Selecting a time window and continuing with the two topoplots was fine.

%%

% this should create 1 figure with 6 subplots
cfg = [];
cfg.layout = layout;
cfg.xlim = linspace(0, 1, 7);
ft_topoplotER(cfg, timelock1);

%%

freq = [];
freq.label = layout.label(1:21); % without COMNT and SCALE
freq.freq = 1:10;
freq.time = (1:50)/50;
freq.powspctrm = randn(21, 10, 50);

%%

cfg = [];
cfg.layout = layout;
ft_multiplotTFR(cfg, freq);

%%

% this should create 1 figure with 6 subplots
cfg = [];
cfg.layout = layout;
cfg.xlim = linspace(0, 1, 7);
cfg.ylim = [4 5];
ft_topoplotTFR(cfg, freq);

%%

comp = [];
comp.topolabel = layout.label(1:21);
comp.label = arrayfun(@num2str, 1:10, 'UniformOutput', false);
comp.topo = randn(21,10);

cfg = [];
cfg.layout = layout;
cfg.component = 1:10;
cfg.zlim = [-2.5 2.5];
ft_topoplotIC(cfg, comp);

