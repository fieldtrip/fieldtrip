function test_bug2620

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_multiplotER ft_multiplotER

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2620.mat'));


close all

cfg = [];
cfg.showoutline = 'yes';
figure; ft_multiplotER(cfg, timelock);

cfg.showoutline = 'yes';
figure; ft_multiplotTFR(cfg, freq);

cfg.showoutline = 'no';
figure; ft_multiplotER(cfg, timelock); % THIS WAS THE PROBLEMATIC ONE

cfg.showoutline = 'no';
figure; ft_multiplotTFR(cfg, freq);
