function test_bug1298

% MEM 3gb
% WALLTIME 00:10:00

% TEST ft_timelockanalysis ft_prepare_leadfield ft_sourceanalysis 

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];
ft_default.feedback = 'no';

megraw = load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/dataFIC.mat'));

cfg = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
megtlock = ft_timelockanalysis(cfg,megraw.dataFIC);

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/vol/Subject01vol_localspheres.mat'))

cfg = [];
cfg.vol = vol;
grid = ft_prepare_leadfield(cfg,megtlock);

cfg = [];
cfg.vol = vol;
cfg.method = 'lcmv';
cfg.grid = grid;
cfg.keepleadfield = 'yes';
cfg.lcmv.keepfilter = 'yes';
megsource1 = ft_sourceanalysis(cfg,megtlock);

cfg = [];
cfg.vol = vol;
cfg.method = 'lcmv';
cfg.grid = grid;
cfg.grid.leadfield = megsource1.leadfield;
cfg.grid.filter = megsource1.avg.filter;
% cfg.keeptrials = 'yes';
cfg.rawtrial = 'yes';
megsource11 = ft_sourceanalysis(cfg,megtlock);

% and then single trial estimates are in megsource11.trial

