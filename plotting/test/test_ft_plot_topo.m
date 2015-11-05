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
isolines = 0:.1:1;

% 1 default style (surfiso) (without isolines is equal to style surf)
figure
ft_plot_topo(x, y, val);
hold on; plot(x, y, 'k*')
axis tight

% 2 style surfiso shading interp
figure
ft_plot_topo(x, y, val, 'style', 'surfiso', 'shading', 'interp');
hold on; plot(x, y, 'k*')
axis tight

% 3 default style (surfiso) with isolines
figure
ft_plot_topo(x, y, val, 'isolines', isolines);
hold on; plot(x, y, 'k*')
axis tight


% 4 style isofill with isolines
figure
ft_plot_topo(x, y, val, 'style', 'isofill', 'isolines', isolines);
hold on; plot(x, y, 'k*')
axis tight



%%% MASKING
% 5 style surf with binary mask
figure
ft_plot_topo(x, y, val, 'datmask', val>0.5, 'style', 'surf', 'clim', [0 1]);
hold on; plot(x, y, 'k*')
axis tight

% 6 style imsat with binary mask
figure
ft_plot_topo(x, y, val, 'datmask', val>0.5, 'style', 'imsat', 'clim', [0 1]);
hold on; plot(x, y, 'k*')
axis tight

% 7 style surf with continuous mask  (part of image will always be masked out because val is always 0 in one spot)
figure
ft_plot_topo(x, y, val, 'datmask', val, 'style', 'surf', 'clim', [0 1]);
hold on; plot(x, y, 'k*')
axis tight

% 8 style imsat with continuous mask  (part of image will always be masked out because val is always 0 in one spot)
figure
ft_plot_topo(x, y, val, 'datmask', val, 'style', 'imsat', 'clim', [0 1]);
hold on; plot(x, y, 'k*')
axis tight

% 9 style imsatiso with continuous mask and isolines (part of image will always be masked out because val is always 0 in one spot)
figure
ft_plot_topo(x, y, val, 'datmask', val, 'style', 'imsatiso', 'clim', [0 1], 'isolines', isolines);
hold on; plot(x, y, 'k*')
axis tight





