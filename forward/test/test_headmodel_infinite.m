% script to test ft_headmodel_infinite. this function is called by
% ft_prepare_headmodel

cfg = [];
cfg.method = 'infinite';
vol1 = ft_prepare_headmodel(cfg);
