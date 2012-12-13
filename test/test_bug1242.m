function test_bug1242

% TEST test_bug1242
% TEST ft_databrowser

load /home/common/matlab/fieldtrip/data/test/bug1242.mat

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, timelock);

