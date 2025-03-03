function test_bug2354

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_topoplotER ft_multiplotER ft_singleplotER
% DATA private

% example ERF data
filename = dccnpath('/project/3031000.02/test/bug2354.mat');
load(filename);

% plot a time window 
cfg = [];
cfg.parameter    = 'individual';
cfg.layout       = 'neuromag306cmb.lay';
cfg.fontsize     = 14;
cfg.zlim = [0 3]*10e-13; 
cfg.xlim = [0.1 0.2];
figure; ft_topoplotER(cfg, data);
figure; ft_singleplotER(cfg, data);
figure; ft_multiplotER(cfg,data);



