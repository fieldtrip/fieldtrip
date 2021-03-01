function test_ft_badchannel

% WALLTIME 00:10:00
% MEM 2gb

dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');

ft_debug off

%%

cfg = [];
cfg.dataset = dataset;
cfg.channel = 'MEG';
cfg.trl = [1 5*3*300 0]; % first 5 trials, each of 3 seconds
cfg.continuous = 'yes';
cfg.demean = 'yes';
data = ft_preprocessing(cfg);

% cfg = [];
% cfg.length = 1;
% data = ft_redefinetrial(cfg, data);

%%

cfg = [];
cfg.method = 'template';
neighbours = ft_prepare_neighbours(cfg, data);

%%

if false
  % this can be used to determine the thresholds
  cfg = [];
  cfg.method = 'summary';
  cfg.neighbours = neighbours;
  data_clean = ft_rejectvisual(cfg, data);
end
%%
% each method has its own optimal threshold

cfg = [];

cfg.method = 'range';
cfg.threshold = 3000e-15;

cfg.method = 'var';
cfg.threshold = 1.6e-25;

cfg.method = 'std';
cfg.threshold = 4.2e-13;

cfg.method = 'max';
cfg.threshold = 1.5e-12;

cfg.method = 'min';
cfg.threshold = -12e-13;

outcfg = ft_badchannel(cfg, data)

cfg = [];
cfg.method = 'weighted';
cfg.neighbours = neighbours;
cfg.senstype = 'meg';
cfg.badchannel = outcfg.badchannel;
data_repaired = ft_channelrepair(cfg, data);

%%
% each method has its own optimal threshold

cfg = [];
cfg.method = 'neighbstdratio';
cfg.threshold = 0.5;

cfg.method = 'neighbcorr';
cfg.threshold = 0.4;

cfg.feedback = 'yes';
cfg.neighbours = neighbours;

cfg.nbdetect = 'any';
outcfg = ft_badchannel(cfg, data)
cfg.nbdetect = 'most';
outcfg = ft_badchannel(cfg, data)
cfg.nbdetect = 'all';
outcfg = ft_badchannel(cfg, data)
cfg.nbdetect = 'median';
outcfg = ft_badchannel(cfg, data)

cfg = [];
cfg.method = 'weighted';
cfg.neighbours = neighbours;
cfg.senstype = 'meg';
cfg.badchannel = outcfg.badchannel;
data_repaired = ft_channelrepair(cfg, data);



