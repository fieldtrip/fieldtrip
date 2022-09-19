function inspect_issue509

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY

%%
% all defaults, should resemble design from Stefan

cfg = [];
data = ft_steadystatesimulation(cfg);

cfg = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, data);

%%
% the following is an attempt to resemble the design from Cassia

cfg = [];
cfg.ntrials = 120;
cfg.level1.condition  = [1 2]; % attend low/high
cfg.level1.gain       = [1.1 1.2];
cfg.level2.condition  = [1 2]; % sequence goes up/down
cfg.level2.gain       = [1.1 1.2];
cfg.level3.condition  = [1 2 3]; % note 1 2 3
cfg.level3.gain       = [1.1 1.2 1.3];
cfg.stimulus1.mode = 'periodic';
cfg.stimulus1.isi = 1/20;
cfg.baseline = 0.1;
cfg.stimulus1.onset = 0;
cfg.stimulus1.onsetjitter = 0;
cfg.stimulus2.mode = 'off';

data = ft_steadystatesimulation(cfg);

cfg = [];
cfg.viewmode = 'vertical';
cfg.channel = [2 3];
ft_databrowser(cfg, data);
