function test_bug2440

% MEM 500mb
% WALLTIME 00:01:00
% test_bug2440

cfg = [];
cfg.method = 'broadband';
cfg.numtrl = 5;
cfg.trllen = 2;

data = ft_freqsimulation(cfg);


cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 8:12;
cfg.t_ftimwin = 4./cfg.foi; 
cfg.toi = 0:0.1:2;
cfg.keeptrials = 'yes';
freq = ft_freqanalysis(cfg, data);

freq.logspctrm = log(freq.powspctrm);

% plot
cfg = [];
cfg.parameter = 'logspctrm';
ft_singleplotTFR(cfg,freq); 

end