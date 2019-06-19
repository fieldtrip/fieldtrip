function test_issue789

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY xdf2fieldtrip ft_read_header ft_read_data ft_read_event

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/xdf/example_EEG_eyeTracking_rigidBody.xdf');

%%
% test the combination of multiple streams, which resamples them to a common time axis

data = xdf2fieldtrip(filename);

cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.ylim = [-20 20];
ft_databrowser(cfg, data);

%%
% test the reading of a single EEG stream

cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);

cfg = [];
cfg.event = ft_read_event(filename); % also read the events
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.ylim = [-20 20];
ft_databrowser(cfg, data);
