function inspect_ft_trialfun_general

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_trialfun_general ft_trialfun_gui ft_trialfun_show

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');

%%

cfg = [];
cfg.dataset = dataset;

%% 266 trials

cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef = [];
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==266);

%% this should show something on screen

cfg.trialfun = 'ft_trialfun_general'; % or ft_trialfun_show
cfg.trialdef = [];
cfg.trialdef.eventtype = '?';
cfgout = ft_definetrial(cfg)

assert(~isfield(cfgout, 'trl'));

%% this should show something on screen

cfg.trialfun = 'ft_trialfun_general'; % or ft_trialfun_show
cfg.trialdef = [];
cfg.trialdef.eventvalue = '?';
cfgout = ft_definetrial(cfg)

assert(~isfield(cfgout, 'trl'));

%% this should show the GUI

cfg.trialfun = 'ft_trialfun_general'; % or ft_trialfun_gui
cfg.trialdef = [];
cfg.trialdef.eventtype = 'gui';
cfgout = ft_definetrial(cfg)

%% this should show the GUI

cfg.trialfun = 'ft_trialfun_general'; % or ft_trialfun_gui
cfg.trialdef = [];
cfg.trialdef.eventvalue = 'gui';
cfgout = ft_definetrial(cfg)

%% this should show the GUI with a subset of events

cfg.trialfun = 'ft_trialfun_general'; % or ft_trialfun_gui
cfg.trialdef = [];
cfg.trialdef.eventtype = 'STIM';
cfg.trialdef.eventvalue = 'gui';
cfgout = ft_definetrial(cfg)

%% one long segment

cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef = [];
cfg.trialdef.length = inf;
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==1);

%% one long segment

cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef = [];
cfg.trialdef.ntrials = 1;
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==1);

%% 266*3 segments

cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef = [];
cfg.trialdef.length = 1;
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==266*3);

%% 266 original trials

cfg.trialfun = 'ft_trialfun_trial';
cfg.trialdef = [];
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==266);

%% ntrials=inf should also work, see https://github.com/fieldtrip/fieldtrip/issues/1912

cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef = [];
cfg.trialdef.length = 1;
cfg.trialdef.ntrials = inf;
cfgout = ft_definetrial(cfg)

assert(size(cfgout.trl,1)==266*3);

