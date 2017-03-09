function test_bug2025

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_plot_vol ft_compute_leadfield

vol = [];
assert(ft_voltype(vol, 'infinite'))

vol.type = 'infinite';

% there is nothing to plot, but nevertheless it should not result in an error
ft_plot_vol(vol);

