function test_issue1507

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY ft_multiplotER ft_singleplotER ft_selectdata
% DATA public


datadir = dccnpath('/project/3031000.02/external/download/tutorial/eventrelatedaveraging');
load(fullfile(datadir, 'dataFC_LP.mat'));
load(fullfile(datadir, 'dataFIC_LP.mat'));

t1 = ft_timelockanalysis([],dataFC_LP);
t2 = ft_timelockanalysis([],dataFIC_LP);

cfg = [];
cfg.latency = [-0.1 1];
t1short = ft_selectdata(cfg, t1);

cfg = [];
cfg.layout = 'CTF151_helmet.mat';

ft_multiplotER(cfg, t2, t1);
ft_multiplotER(cfg, t2, t1short);

cfg.select = 'union';
ft_multiplotER(cfg, t2, t1short);


cfg = rmfield(cfg, 'select');
cfg.channel = 'MLT';
ft_singleplotER(cfg, t2, t1);
ft_singleplotER(cfg, t2, t1short);
cfg.select = 'union';
ft_singleplotER(cfg, t2, t1short);
ft_singleplotER(cfg, t2);
