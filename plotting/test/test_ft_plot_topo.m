function test_ft_plot_topo

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_topo
% TEST ft_plot_topo

x = randn(30,1);
y = randn(30,1);

% make something round between 0 and 1
val = -sqrt(x.^2 + y.^2);
val = val-min(val);
val = val./max(val);

figure
ft_plot_topo(x, y, val);
hold on; plot(x, y, 'k*')

figure
ft_plot_topo(x, y, val, 'shading', 'interp');
hold on; plot(x, y, 'k*')

figure
ft_plot_topo(x, y, val, 'datmask', val>0.5, 'style', 'surf');
hold on; plot(x, y, 'k*')

figure
ft_plot_topo(x, y, val, 'datmask', val>0.5, 'style', 'imsat', 'clim', [0 1]);
hold on; plot(x, y, 'k*')



