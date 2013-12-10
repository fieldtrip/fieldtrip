function test_ft_plot_box

% MEM 1gb
% WALLTIME 0:03:00

% TEST test_ft_plot_box
% TEST ft_plot_box

figure
ft_plot_box([-1 1 2 3], 'facecolor', 'b');
axis([-4 4 -4 4])
