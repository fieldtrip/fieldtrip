function test_ft_sourcestatistics

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_sourcestatistics

gridsize = 1320;
nsamples = 20;
source = [];
source.pos = randn(gridsize,3);
source.time = 1:nsamples;
source.pow = randn(gridsize,nsamples);
source.inside = 1:gridsize/2;

source2 = source;
source2.pow = randn(gridsize,nsamples);

source3 = source;
source3.pow = randn(gridsize,nsamples);

% test case of Monte Carlo method
cfg = [];
cfg.parameter = 'pow';
cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.numrandomization = 'all';
cfg.design = [1 -0.5 -0.5];
sourceout = ft_sourcestatistics(cfg, source, source2, source3);

