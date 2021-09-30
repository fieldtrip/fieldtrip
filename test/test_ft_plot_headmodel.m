function test_ft_plot_headmodel

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_plot_headmodel

vol.r = [86 88 92 100];
vol.o = [0 0 40];

figure; ft_plot_headmodel(vol, 'facecolor', 'skin_light');
figure; ft_plot_headmodel(vol, 'facecolor', 'skin_medium_light');
figure; ft_plot_headmodel(vol, 'facecolor', 'skin_medium');
figure; ft_plot_headmodel(vol, 'facecolor', 'skin_medium_dark');
figure; ft_plot_headmodel(vol, 'facecolor', 'skin_dark');
