function test_ft_timelockstatistics

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_timelockstatistics, findcluster, clusterstat, ft_statistics_montecarlo

%For the case of "chan_time"

% make fake dataset
timelock = cell(1,10);
for idat = 1:10
  timelock{idat}.label = {'chan1','chan2','chan3'};
  timelock{idat}.dimord = 'chan_time';
  timelock{idat}.avg = rand(3,30);
  timelock{idat}.time = 0.1:0.1:3;
  timelock{idat}.cfg = [];
end

% do stats - montecarlo
cfg = [];
neighbours(1).label = 'chan1';
neighbours(1).neighblabel = {'chan2', 'chan3'};
neighbours(2).label = 'chan2';
neighbours(2).neighblabel = {'chan1', 'chan3'};
neighbours(3).label = 'chan3';
neighbours(3).neighblabel = {'chan1', 'chan2'};
cfg.neighbours  = neighbours;
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05; 
cfg.correctm    = 'cluster'; 
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'parametric';
cfg.numrandomization = 500;
cfg.design = [ones(1,5) ones(1,5).*2; 1:5 1:5;];
cfg.ivar   = 1;
cfg.uvar   = 2;
stat = ft_timelockstatistics(cfg,timelock{:});

% do stats - analytic
cfg = [];
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05; 
cfg.design = [ones(1,5) ones(1,5).*2; 1:5 1:5;];
cfg.ivar   = 1;
cfg.uvar   = 2;
stat = ft_timelockstatistics(cfg,timelock{:});


% do stats - analytic
cfg = [];
cfg.method      = 'stats';
cfg.statistic   = 'ttest';
cfg.alpha       = 0.05; 
cfg.design = [ones(1,10) ];
stat = ft_timelockstatistics(cfg,timelock{:});


