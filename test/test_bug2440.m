function test_bug2440

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_bug2440 ft_freqsimulation ft_freqanalysis ft_singleplotTFR

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
cfg.trials = 2;
ft_singleplotTFR(cfg,freq); 

end
