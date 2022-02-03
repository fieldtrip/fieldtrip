load('some_trial_data.mat')

cfg = [];
cfg.artfctdef.eog.channel = 'VEOG';
cfg.artfctdef.eog.cutoff = 2;
% cfg.artfctdef.eog.interactive = 'yes';
cfg.continuous = 'no';
cfg = ft_artifact_eog(cfg,data);
blinks = cfg.artfctdef.eog;

cfg = [];
cfg.event = event;
cfg.ylim = [-100 100];
cfg.artfctdef.eog = blinks;
cfg.artfctdef.visual.artifact = [];

cfg2 = ft_databrowser_multitrial(cfg,data);

cfg3 = ft_databrowser(cfg2,data);

%%
ncfg = [];
ncfg.event = event;
ncfg.ylim = [-100 100];
ncfg.artfctdef.eog = blinks;
ncfg.blocksize = 8.504;
ncfg.continuous = 'noI';
ncfg.ploteventlabels = 'no';

cfg4 = ft_databrowser(ncfg,data);
