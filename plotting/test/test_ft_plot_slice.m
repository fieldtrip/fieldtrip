function test_ft_plot_slice

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_slice
% TEST ft_plot_slice

dat = randn(64, 64, 64);

figure
ft_plot_slice(dat);

