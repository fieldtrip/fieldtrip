function test_pull2099

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY ft_multiplotER ft_singleplotER ft_topoplotER lineattributes_common ft_databrowser
% PRIVATEDATA


%%
load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull2099.mat'));

% try a bunch of visualizations
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
ft_multiplotER(cfg, t);

cfg = [];
cfg.linecolor = 'spatial';
cfg.layout = 'butterfly';
ft_multiplotER(cfg, t); % ->monochrome because spatial information of the sensors is not known (grad missing)

cfg = [];
cfg.linecolor = 'spatial';
cfg.layout = 'CTF275_helmet.mat';
cfg.linewidth = 2;
ft_multiplotER(cfg, t); 

cfg = [];
cfg.linecolor = 'spatial';
cfg.viewmode = 'butterfly';
cfg.layout = 'CTF275_helmet.mat';
ft_multiplotER(cfg, t);

cfg = [];
cfg.viewmode = 'butterfly';
cfg.layout = 'CTF275_helmet.mat';
cfg.colorgroups = 'labelchar3';
ft_multiplotER(cfg, t);

cfg = [];
cfg.channel = 'MRO';
cfg.showlocations = 'yes';
cfg.linecolor = [0.8 0.4 0.3];
cfg.layout = 'CTF275_helmet.mat';
ft_singleplotER(cfg, t);

cfg = [];
cfg.channel = 'MRO';
cfg.showlocations = 'yes';
cfg.layout = 'CTF275_helmet.mat';
ft_singleplotER(cfg, t1, t2, t3, t4);



