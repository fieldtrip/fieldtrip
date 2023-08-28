function test_bug2476

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqanalysis ft_freqanalysis_mvar
% DATA private

cfg           = [];
cfg.inputfile = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2476.mat');
cfg.method    = 'mvar';
freq          = ft_freqanalysis(cfg);
