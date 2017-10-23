function test_bug2476

% WALLTIME 00:10:00
% MEM 150mb

% TEST ft_freqanalysis ft_freqanalysis_mvar

cfg           = [];
cfg.inputfile = dccnpath('/home/common/matlab/fieldtrip/data/test/test_bug2476.mat');
cfg.method    = 'mvar';
freq          = ft_freqanalysis(cfg);
