function test_issue1217

% MEM 4gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivityanalysis ft_laggedcoherence fourierspctrm2lcrsspctrm

% see https://github.com/fieldtrip/fieldtrip/issues/1217
% and https://github.com/fieldtrip/fieldtrip/pull/1233

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1217.mat'));

% do a spectral decomposition first
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'fourier';
cfg.foi    = 2:2:100;
cfg.t_ftimwin = ones(1,numel(cfg.foi))./2;
cfg.taper  = 'hanning';
cfg.toi    = -0.5:0.05:1.5;
freq = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'laggedcoherence';
cfg.laggedcoherence.lags = 0.5;
lcoh = ft_connectivityanalysis(cfg, freq);

% try it on a single frequency
cfg = [];
cfg.frequency = 10;
freq2 = ft_selectdata(cfg, freq);

cfg = [];
cfg.method = 'laggedcoherence';
cfg.laggedcoherence.lags = 0.5;
lcoh2 = ft_connectivityanalysis(cfg, freq2);
