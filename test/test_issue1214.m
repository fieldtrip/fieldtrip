function test_issue1214

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_singleplotER ft_multiplotER ft_plot_vector

%% make some data

timelock1 = ft_timelockanalysis([], ft_freqsimulation([]));

cfg = [];
cfg.parameter = 'avg';
cfg.operation = '-1*x1';
timelock2 = ft_math(cfg, timelock1);

%% this works

cfg = [];
cfg.channel = 'mix';
figure; ft_singleplotER(cfg, timelock1, timelock2);

%% this works

cfg = [];
cfg.channel = 'mix';
cfg.graphcolor = 'my';
figure; ft_singleplotER(cfg, timelock1, timelock2);

%% this failed

cfg = [];
cfg.channel = 'mix';
cfg.graphcolor = [
  1 0 0
  0 1 0
  ];
figure; ft_singleplotER(cfg, timelock1, timelock2);

%% this failed

cfg = [];
cfg.layout = 'vertical';
cfg.graphcolor = [
  1 0 0
  0 1 0
  ];
figure; ft_multiplotER(cfg, timelock1, timelock2);