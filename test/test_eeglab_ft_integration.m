function test_eeglab_ft_integration

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_eeglab_ft_integration
% TEST ft_dipolefitting ft_checkdata

% See http://bugzilla.fcdonders.nl/show_bug.cgi?id=2595

load test_EEGLAB_ft_integration.mat

%% test the original reported problem
source1 = ft_dipolefitting(cfg, data);

%% also do the nonlinear fit
cfg.component = [1 2 3];
cfg.nonlinear = 'yes';
source2 = ft_dipolefitting(cfg, data);

% do not use the optimalization toolbox (which is the default if available)
cfg.dipfit.optimfun = 'fminsearch';
source3 = ft_dipolefitting(cfg, data);
