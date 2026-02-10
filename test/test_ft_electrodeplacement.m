function test_ft_electrodeplacement

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_electrodeplacement
% DATA no

[headshape.pos, headshape.tri] = mesh_sphere(1000);
headshape.pos = 100*headshape.pos;
headshape.unit = 'mm';
headshape.coordsys = 'ctf';

figure
ft_plot_mesh(headshape)
ft_headlight

%%

nas = [100 0 0];
lpa = [0 100 0];
rpa = [0 -100 0];
ini = [-100 0 0];

cfg = [];
cfg.method = '1020';
cfg.fiducial.nas = nas;
cfg.fiducial.ini = ini;
cfg.fiducial.lpa = lpa;
cfg.fiducial.rpa = rpa;
elec1020 = ft_electrodeplacement(cfg, headshape);

figure
ft_plot_mesh(headshape, 'facecolor', 'skin');
ft_plot_sens(elec1020, 'headshape', headshape, 'elecshape', 'disc');
ft_headlight

%%

cfg = [];
cfg.method = 'equidistant';
cfg.numelec = 60;
cfg.fiducial.nas = [100 0 0];
cfg.fiducial.ini = [-100 0 0];
cfg.fiducial.lpa = [0 100 0];
cfg.fiducial.rpa = [0 -100 0];
cfg.feedback = 'yes';
equidistant = ft_electrodeplacement(cfg, headshape);

figure
ft_plot_mesh(headshape, 'facecolor', 'skin', 'edgecolor', 'none', 'axes', 'on'); alpha 0.5
ft_plot_sens(equidistant, 'headshape', [], 'label', 'label', 'elecshape', 'sphere');
ft_headlight
