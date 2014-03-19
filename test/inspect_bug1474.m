function inspect_bug1474

% TEST inspect_bug1474
% TEST ft_rejectvisual

% this script tests ft_rejectvisual with cfg.method='summary' for the case
% where the data only contains one channel

% since ft_rejectvisual requires user interaction with the GUI, we should
% not run this automatically

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1474.mat

cfg = [];
cfg.method = 'summary';
ft_rejectvisual(cfg, data);
