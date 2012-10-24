function test_bug1474

% TEST test_bug1474
% TEST ft_rejectvisual

% this script tests ft_rejectvisual with cfg.method='summary' for the case
% where the data only contains one channel

% since ft_rejectvisual requires user interaction with the GUI, we should
% just return here, to avoid false negatives in the testing suite.
% The rest of the code is only for reference.
return;

load test_bug1474.mat;

cfg = [];
cfg.method = 'summary';
ft_rejectvisual(cfg, data);