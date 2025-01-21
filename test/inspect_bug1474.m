function inspect_bug1474

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_rejectvisual
% DATA private

% this script tests ft_rejectvisual with cfg.method='summary' for the case
% where the data only contains one channel

% since ft_rejectvisual requires user interaction with the GUI, we should
% not run this automatically

cd(dccnpath('/project/3031000.02/test'))
load bug1474.mat

cfg = [];
cfg.method = 'summary';
ft_rejectvisual(cfg, data);
