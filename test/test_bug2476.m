function test_bug2476

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqanalysis ft_freqanalysis_mvar

cfg           = [];
cfg.inputfile = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2476.mat');
cfg.method    = 'mvar';
freq          = ft_freqanalysis(cfg);
