function test_headmodel_infinite

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_headmodel ft_headmodel_infinite

% function to test ft_headmodel_infinite. this function is called by
% ft_prepare_headmodel

cfg = [];
cfg.method = 'infinite';
vol1 = ft_prepare_headmodel(cfg);
