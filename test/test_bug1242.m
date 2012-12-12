function test_bug1242

% TEST test_bug1242
% TEST ft_databrowser

load test_bug1242;

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, timelock);

