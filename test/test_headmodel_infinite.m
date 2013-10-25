function test_headmodel_infinite

% MEM 1gb
% WALLTIME 00:03:01

% TEST test_headmodel_infinite
% TEST ft_prepare_headmodel ft_headmodel_infinite

% function to test ft_headmodel_infinite. this function is called by
% ft_prepare_headmodel

cfg = [];
cfg.method = 'infinite';
vol1 = ft_prepare_headmodel(cfg);
