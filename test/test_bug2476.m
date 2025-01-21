function test_bug2476

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqanalysis ft_freqanalysis_mvar
% DATA private

cfg           = [];
cfg.inputfile = dccnpath('/project/3031000.02/test/bug2476.mat');
cfg.method    = 'mvar';
freq          = ft_freqanalysis(cfg);
