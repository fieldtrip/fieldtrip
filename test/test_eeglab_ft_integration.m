function test_eeglab_ft_integration

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_eeglab_ft_integration
% TEST ft_dipolefitting ft_checkdata

% See http://bugzilla.fcdonders.nl/show_bug.cgi?id=2595

load test_EEGLAB_ft_integration.mat

source = ft_dipolefitting(cfg, data);
