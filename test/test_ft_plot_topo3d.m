function test_ft_plot_topo3d

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_topo3d
% TEST ft_plot_topo3d

pos = randn(10,3);
val = randn(10,1);

figure
ft_plot_topo3d(pos, val);
