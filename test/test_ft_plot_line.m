function test_ft_plot_line

% MEM 1gb
% WALLTIME 00:10:00

% DEPENDENCY ft_plot_line

x = randn(30,1);
y = randn(30,1);

figure
ft_plot_line(x, y);
