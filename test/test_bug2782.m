function test_bug2782

% WALLTIME 00:10:00
% MEM 1gb

timelock = [];
timelock.label = {'1', '2', '3'};
timelock.time = (1:1000)/1000;
timelock.avg = randn(3,1000);
timelock.trial = randn(1,3,1000);
timelock.dimord = 'rpt_chan_time';
timelock.trialinfo = 1:11;

cfg = [];
cfg.trials = 1;
tmp = ft_selectdata(cfg,timelock)

% this should not change
assert(size(tmp.trialinfo,1)==1);
assert(size(tmp.trialinfo,1)==11);



