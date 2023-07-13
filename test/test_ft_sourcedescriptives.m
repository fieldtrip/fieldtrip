function test_ft_sourcedescriptives

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_sourcedescriptives
% DATA no

gridsize = 38*48*41;
nsamples = 20;
source = [];
source.pos = randn(gridsize,3);
source.time = 1:nsamples;
source.pow = randn(gridsize,nsamples);
source.inside = 1:gridsize;

cfg = [];
cfg.powmethod = 'regular';
sourceout = ft_sourcedescriptives(cfg,source);
