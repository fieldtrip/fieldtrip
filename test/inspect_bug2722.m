function inspect_bug2722

% WALLTIME 00:10:00
% MEM 1gb

% TEST test_bug2721
% TEST ft_multiplotTFR

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2722.mat'));

%%

cfg = [];
cfg.layout        = 'neuromag306mag.lay';

figure
ft_multiplotER(cfg, timelock_dif);

%%

timelock_std.mask = timelock_dif.avg>2e-13; % this one is used in multiplot and singleplot
timelock_dev.mask = timelock_dif.avg>2e-13; % both are used in topoplot

cfg.maskparameter = 'mask';
cfg.maskstyle     = 'box';

figure
ft_multiplotER(cfg, timelock_std, timelock_dev);

% INSTRUCTION: use the interactive mode to make a selection in the multiplot and
% the singleplot
