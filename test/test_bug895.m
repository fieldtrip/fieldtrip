function test_bug895

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY

% using statfun_indepsamplesZcoh results in an output structure that
% contains 'chan' in the dimord, rather than 'chancmb'

% load some data
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275.mat'));

% do spectral transformation
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper  = 'hanning';
cfg.foilim = [1 20];
freq = ft_freqanalysis(cfg, data);

% compute statistic
cfg = [];
cfg.method    = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesZcoh';
cfg.parameter = 'fourierspctrm';
cfg.numrandomization = 1;
cfg.design    = [ones(1,5) ones(1,5)*2];
stat = ft_freqstatistics(cfg, freq);

assert(strcmp(stat.dimord,'chancmb_freq'));
