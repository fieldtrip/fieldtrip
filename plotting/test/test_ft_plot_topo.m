function test_ft_plot_topo

% MEM 1gb
% WALLTIME 0:03:00

% TEST test_ft_plot_topo
% TEST ft_plot_topo

x = randn(10,1);
y = randn(10,1);
val = randn(10,1);

figure
ft_plot_topo(x, y, val);
