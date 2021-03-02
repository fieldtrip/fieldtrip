function test_ft_badchannel

% WALLTIME 00:10:00
% MEM 2gb

dataset = dccnpath('/home/common/matlab/fieldtrip/datatrl/Subject01.ds');

ft_debug off

%%

cfg = [];
cfg.dataset = dataset;
cfg.channel = 'MEG';
cfg.trl = [1 7*3*300 0]; % first 7 trials, each of 3 seconds
cfg.continuous = 'yes';
cfg.demean = 'yes';
datacnt = ft_preprocessing(cfg);

cfg = [];
cfg.length = 1;
datatrl = ft_redefinetrial(cfg, datacnt);

%%

cfg = [];
cfg.method = 'template';
neighbours = ft_prepare_neighbours(cfg, datatrl);

%%

% this can be used to determine the thresholds
cfg = [];
cfg.method = 'summary';
cfg.neighbours = neighbours;
data_clean = ft_rejectvisual(cfg, datatrl);

%%
% each method has its own optimal threshold
% FT_REJECTVISUAL can be used to determine the appropriate value



cfg = [];
cfg.feedback = 'yes';
cfg.neighbours = neighbours;

% cfg.method = 'range';
% cfg.threshold = 3000e-15;
%
% cfg.method = 'var';
% cfg.threshold = 1.6e-25;
%
% cfg.method = 'std';
% cfg.threshold = 4.2e-13;
%
% cfg.method = 'max';
% cfg.threshold = 1.5e-12;
%
% cfg.method = 'min';
% cfg.threshold = -12e-13;
%
% cfg.method = 'min';
% cfg.threshold = -12e-13;
%
% cfg.method = 'zvalue';
% cfg.threshold = 3;

cfg.method = 'neighbexpvar';
cfg.threshold = 0.5;

% cfg.method = 'neighbstdratio';
% cfg.threshold = 0.5;
%
% cfg.method = 'neighbcorr';
% cfg.threshold = 0.5;
%
% these options can be varied for neighbstdratio and neighbcorr
% cfg.nbdetect = 'median';
% cfg.nbdetect = 'any';
% cfg.nbdetect = 'most';
% cfg.nbdetect = 'all';

outcfg = ft_badchannel(cfg, datatrl)

%%

cfg = [];
cfg.method = 'weighted';
cfg.neighbours = struct([]); % neighbours;
cfg.senstype = 'meg';
cfg.badchannel = outcfg.badchannel;
data_repaired = ft_channelrepair(cfg, datatrl);


