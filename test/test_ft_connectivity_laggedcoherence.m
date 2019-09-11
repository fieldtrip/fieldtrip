function test_ft_connectivity_laggedcoherence

% MEM 1500mb
% WALLTIME 00:10:00

% DEPENDENCY: test_ft_connectivity_laggedcoherence
% DEPENDENCY: ft_connectivityanalysis
% DEPENDENCY: fourierspctrm2lcrsspctrm

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1217.mat'));

% do a spectral decomposition first
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'fourier';
cfg.foi    = 2:2:100
cfg.t_ftimwin = ones(1,numel(cfg.foi))./2;
cfg.taper  = 'hanning';
cfg.toi    = -0.5:0.05:1.5;
freq = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'laggedcoherence';
cfg.laggedcoherence.lags = 0.5;
lcoh = ft_connectivityanalysis(cfg, freq);
