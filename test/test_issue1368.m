function test_issue1368

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_plot_sens mesh_cylinder

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1368/hull.mat'))
load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1368/elec.mat'))

%%

ft_plot_sens(elec, 'elecshape', 'point', 'headshape', []);
ft_plot_sens(elec, 'elecshape', 'point', 'headshape', hull);

%%

ft_plot_sens(elec, 'elecshape', 'circle', 'headshape', []);
ft_plot_sens(elec, 'elecshape', 'circle', 'headshape', hull);

%%

ft_plot_sens(elec, 'elecshape', 'square', 'headshape', []);
ft_plot_sens(elec, 'elecshape', 'square', 'headshape', hull);

%%

ft_plot_sens(elec, 'elecshape', 'sphere', 'headshape', []);
ft_plot_sens(elec, 'elecshape', 'sphere', 'headshape', hull);

%%

ft_plot_sens(elec, 'elecshape', 'disc', 'headshape', []);
ft_plot_sens(elec, 'elecshape', 'disc', 'headshape', hull);
