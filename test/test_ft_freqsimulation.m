function test_ft_freqsimulation

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_freqsimulation

cfg = [];
cfg.method = 'superimposed';
dataout = ft_freqsimulation(cfg);