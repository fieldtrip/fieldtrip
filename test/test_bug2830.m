function test_bug2830

% WALLTIME 00:20:00
% MEM 3gb

% TEST ft_sourcestatistics ft_statistics_montecarlo clusterstat

% the test directory holds the data and the statfun
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2830'))

load data_bug2830.mat

%% this is the original code causing the problem

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

%% pretend it is organized as a slab from a 3-D volume

for i=1:44
  allsource{i}.dim = [22*27 23 1]; % note, still 3-D
end  

statpow = ft_sourcestatistics(cfg, allsource{:,1}, allsource{:,2});

%% pretend it is organized as a strip of voxels from a 3-D volume

for i=1:44
  allsource{i}.dim = [22*27*23 1 1]; % note, still 3-D
end  

statpow = ft_sourcestatistics(cfg, allsource{:,1}, allsource{:,2});

%% add a time and frequency dimension with multiple elements

for i=1:44
  allsource{i}.dim = [22 27 23];
  allsource{i}.freq = [30 40 50 60];
  allsource{i}.time = [1 2];
  allsource{i}.avg.pow = randn([22*27*23 4 2]); % pos_freq_time
end  

cfg.numrandomization = 100; % otherwise it takes so long
statpow = ft_sourcestatistics(cfg, allsource{:,1}, allsource{:,2});

