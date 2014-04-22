function test_ft_plot_vol

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_vol
% TEST ft_plot_vol

vol.r = [86 88 92 100];
vol.o = [0 0 40];
figure
ft_plot_vol(vol);
