% function test_bug2773

% WALLTIME 00:20:00
% MEM 2gb

% TEST test_bug2773
% TEST ft_dipolefitting ft_movieplotER ft_prepare_sourcemodel


orig = load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2773.mat'));
vol  = orig.cfg.vol;
elec = orig.cfg.elec;

cfg = [];
cfg.elec = elec;
layout = ft_prepare_layout(cfg);

figure; ft_plot_layout(layout);
figure; ft_plot_sens(elec);
figure; ft_plot_vol(vol);




