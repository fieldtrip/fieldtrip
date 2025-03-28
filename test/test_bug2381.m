function test_bug2381

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_sourcestatistics
% DATA private

%filename = dccnpath(fullfile('/project/3031000.02/test/bug2381','AVG.mat'));
filename = dccnpath('/project/3031000.02/test/bug2381.mat');
load(filename);

nsubj = 15;
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.parameter = 'pow';
cfg.correctm  = 'cluster';
cfg.numrandomization = 10;
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj) 2*ones(1,nsubj)];
cfg.uvar = 1;
cfg.ivar = 2;
%stat = ft_sourcestatistics(cfg, AverageSource1, AverageSource2);
stat = ft_sourcestatistics(cfg, source1, source2);
