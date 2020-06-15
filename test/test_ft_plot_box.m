function test_ft_plot_box

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_plot_box

figure
ft_plot_box([-1 1 2 3], 'facecolor', 'b');
axis([-4 4 -4 4])
