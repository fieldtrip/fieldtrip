function test_ft_plot_matrix

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_matrix
% TEST ft_plot_matrix

dat = randn(30,50);

figure
ft_plot_matrix(dat);
