function test_laggedcoherence

% MEM 1500mb
% WALLTIME 00:20:00
% DEPENDENCY ft_connectivityanalysis ft_laggedcoherence ft_freqanalysis

% this is a test for the low-level function

datain = [];
datain.label    = {'Ch1','Ch2'};
datain.fsample  = 1000;
datain.time{1}  = 0:1/datain.fsample:1.832;
datain.trial{1} = randn(2,length(datain.time{1}));

cfg=[];
cfg.foi = [5 6 7 8 17 18 19 20 21];
cfg.loi = [1 2 3];
cfg.numcycles = 3;
lcout = ft_laggedcoherence(cfg,datain);

