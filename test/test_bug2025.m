function test_bug2025

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_plot_headmodel ft_compute_leadfield

vol = [];
assert(ft_headmodeltype(vol, 'infinite'))

vol.type = 'infinite';

% there is nothing to plot, but nevertheless it should not result in an error
ft_plot_headmodel(vol);

