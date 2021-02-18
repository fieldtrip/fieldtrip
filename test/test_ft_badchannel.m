function test_ft_badchannel

% WALLTIME 00:10:00
% MEM 2gb

dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');

%%

cfg = [];
cfg.dataset = dataset;
cfg.channel = 'MEG';
cfg.trl = [1 5*3*300 0]; % first 5 trials, each of 3 seconds
cfg.continuous = 'yes';
data = ft_preprocessing(cfg);

%%

cfg = [];
cfg.method = 'template';
neighbours = ft_prepare_neighbours(cfg, data);

%%

cfg = [];
cfg.method = 'std';
cfg.threshold = 400e-15;
outcfg = ft_badchannel(cfg, data)

cfg = [];
cfg.method = 'weighted';
cfg.neighbours = neighbours;
cfg.senstype = 'meg';
cfg.badchannel = outcfg.badchannel;
data_repaired = ft_channelrepair(cfg, data);

%%

cfg = [];
cfg.method = 'range';
cfg.threshold = 3000e-15;
outcfg = ft_badchannel(cfg, data)

cfg = [];
cfg.method = 'weighted';
cfg.neighbours = neighbours;
cfg.senstype = 'meg';
cfg.badchannel = outcfg.badchannel;
data_repaired = ft_channelrepair(cfg, data);

%%

cfg = [];
cfg.method = 'neighstdratio';
cfg.threshold = 0.5;
cfg.neighbours = neighbours;
outcfg = ft_badchannel(cfg, data)

cfg = [];
cfg.method = 'weighted';
cfg.neighbours = neighbours;
cfg.senstype = 'meg';
cfg.badchannel = outcfg.badchannel;
data_repaired = ft_channelrepair(cfg, data);

%%

cfg = [];
cfg.method = 'neighcorrel';
cfg.threshold = 0.4;
cfg.neighbours = neighbours;
outcfg = ft_badchannel(cfg, data)

cfg = [];
cfg.method = 'weighted';
cfg.neighbours = neighbours;
cfg.senstype = 'meg';
cfg.badchannel = outcfg.badchannel;
data_repaired = ft_channelrepair(cfg, data);


