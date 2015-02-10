function test_bug2830

% WALLTIME = 00:10:00
% MEM = 2gb

% TEST test_bug2830
% TEST ft_sourcestatistics ft_statistics_montecarlo clusterstat

% the test directory holds the data and the statfun
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2830'))

load data_bug2830.mat

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'freq_statfun_relchange';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05 ;
cfg.clusterstatistic = 'max';
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 1000;

cfg.clusterthreshold = 'nonparametric_individual';

cfg.design = [
  ones(1,size(allsource,1)) ones(1,size(allsource,1))*2;
  1:size(allsource,1) 1:size(allsource,1)
  ];
cfg.ivar = 1;
cfg.uvar = 2;

cfg.parameter = 'avg.pow';

statpow = ft_sourcestatistics(cfg, allsource{:,1}, allsource{:,2});
