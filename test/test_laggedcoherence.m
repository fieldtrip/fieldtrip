function test_laggedcoherence

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_connectivityanalysis ft_laggedcoherence ft_freqanalysis

% this is a test for the low-level function

datain=[];
datain.label={'Ch1','Ch2'};
datain.fsample=1000;
datain.time{1}=0:1/datain.fsample:8.832; % changed by JM, since it does not make sense to have less than 6 s of data if 1 Hz is requested with 3 cycles
datain.trial{1}=randn(2,length(datain.time{1}));

cfg=[];
cfg.foi=[1 1.5 2 3 4 5 6 7 8 17 18 19 20 21];
cfg.loi=[1 2 3];
cfg.numcycles=3;

lcout=ft_laggedcoherence(cfg,datain);
