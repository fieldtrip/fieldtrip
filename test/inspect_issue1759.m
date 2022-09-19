function inspect_issue1759

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY resampledesign

% resampledesign should give warning if the number of user-selected permutations is not optimal
% this is actually a duplicate of https://github.com/fieldtrip/fieldtrip/issues/1313

timelock = [];
timelock.label = {'1', '2', '3'};
timelock.time = (0:999)/1000 - 0.1;
timelock.dimord = 'chan_time';

nsubj = 10; % there are 2^10 = 1024 possible permutations
condition1 = cell(1,nsubj);
condition2 = cell(1,nsubj);

for i=1:nsubj
  condition1{i} = timelock;
  condition1{i}.avg = randn(3, 1000);
  condition2{i} = timelock;
  condition2{i}.avg = randn(3, 1000);
  condition2{i}.avg(1,301:400) = randn(1, 100) + 0.1;
end

neighbours = [];
neighbours(1).label = '1';
neighbours(1).neighblabel = '2';
neighbours(2).label = '2';
neighbours(2).neighblabel = '3';
neighbours(3).label = '3';
neighbours(3).neighblabel = '1';

%%

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.correcttail = 'alpha';
cfg.tail = 0;
cfg.alpha = 0.025; % since we are doing two tests
cfg.neighbours = neighbours;
cfg.design = [
  1:nsubj         1:nsubj
  1*ones(1,nsubj) 2*ones(1,nsubj)
  ];
cfg.uvar = 1;
cfg.ivar = 2;

%%
% this should give a warning
cfg.numrandomization = 2000;
stat = ft_timelockstatistics(cfg, condition1{:}, condition2{:});

%%
% this should give a warning
cfg.numrandomization = 1000;
stat = ft_timelockstatistics(cfg, condition1{:}, condition2{:});

%%
% this should not give a warning
cfg.numrandomization = 'all';
stat = ft_timelockstatistics(cfg, condition1{:}, condition2{:});

%%
% this should not give a warning
cfg.numrandomization = 500;
stat = ft_timelockstatistics(cfg, condition1{:}, condition2{:});
