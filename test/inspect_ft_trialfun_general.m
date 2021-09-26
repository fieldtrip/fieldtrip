function inspect_ft_trialfun_general

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_trialfun_general ft_trialfun_gui ft_trialfun_show

dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');

%%

cfg = [];
cfg.dataset = dataset;
cfg.trialfun = 'ft_trialfun_general';

%% this should show something on screen

cfg.trialdef = [];
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==266);

%% this should show something on screen

cfg.trialdef = [];
cfg.trialdef.eventtype = 'show';
cfgout = ft_definetrial(cfg)

assert(~isfield(cfgout, 'trl'));

%% this should show something on screen

cfg.trialdef = [];
cfg.trialdef.eventvalue = 'show';
cfgout = ft_definetrial(cfg)

assert(~isfield(cfgout, 'trl'));

%% this should show the GUI

cfg.trialdef = [];
cfg.trialdef.eventtype = 'gui';
cfgout = ft_definetrial(cfg)

%% this should show the GUI

cfg.trialdef = [];
cfg.trialdef.eventvalue = 'gui';
cfgout = ft_definetrial(cfg)

%% this should show the GUI with a subset of events

cfg.trialdef = [];
cfg.trialdef.eventtype = 'STIM';
cfg.trialdef.eventvalue = 'gui';
cfgout = ft_definetrial(cfg)

%% one long segment

cfg.trialdef = [];
cfg.trialdef.length = inf;
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==1);

%% one long segment

cfg.trialdef = [];
cfg.trialdef.ntrials = 1;
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==1);

%% 266*3 segments

cfg.trialdef = [];
cfg.trialdef.length = 1;
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==266*3);

