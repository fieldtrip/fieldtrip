function test_bug2620

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_multiplotER ft_multiplotER
% DATA private

load(dccnpath('/project/3031000.02/test/bug2620.mat'));


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
