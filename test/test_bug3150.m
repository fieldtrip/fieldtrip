function test_bug3150

% WALLTIME 00:10:00
% MEM 2gb

condition1 = {};
condition2 = {};
nsamples = 100000;

%% create the data

for i=1:14
  timelock = [];
  timelock.dimord = 'chan_time';
  timelock.label = {'1'};
  timelock.time = 1:nsamples;
  timelock.avg = 1+randn(1,nsamples);
  
  condition1{i} = timelock;
end

for i=1:19
  timelock = [];
  timelock.dimord = 'chan_time';
  timelock.label = {'1'};
  timelock.time = 1:nsamples;
  timelock.avg = 1+randn(1,nsamples);
  
  condition2{i} = timelock;
end

%% do the statistics
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'indepsamplesT';
cfg.correctm = 'no';
cfg.design = [1*ones(1,14) 2*ones(1,19)];
cfg.ivar = 1;
stat = ft_timelockstatistics(cfg, condition1{:}, condition2{:});

%% this is not part of the bug report, but easy to check here
assert(mean(stat.mask)>0.05-0.02, 'not enough false alarms');
assert(mean(stat.mask)<0.05+0.02, 'too many false alarms');
