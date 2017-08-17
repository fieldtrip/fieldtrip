function inspect_bug1830

% WALLTIME 00:10:00
% MEM 1gb

%%

[v, p] = ft_version;

skin = ft_read_headshape(fullfile(p, 'template', 'headmodel', 'skin', 'standard_skin_14038.vol'));
elec = ft_read_sens(fullfile(p, 'template', 'electrode', 'standard_1020.elc'));

skin.coordsys = 'mni';

% r = rotate([0 0 15])
r = [
  0.9659   -0.2588         0         0
  0.2588    0.9659         0         0
  0         0    1.0000         0
  0         0         0    1.0000
  ];


% t = translate([10 20 30]);
t = [
  1     0     0    10
  0     1     0    10
  0     0     1    20
  0     0     0     1
  ];

elec = ft_transform_sens(t * r, elec);

% figure
% ft_plot_mesh(skin, 'facecolor', 'skin');
% ft_plot_sens(elec);
% alpha 0.5


%%

cfg = [];
cfg.template.headshape = skin;
cfg.template.headshapestyle = 'surface';

cfg.individual.elec = elec;
out = ft_interactiverealign(cfg);

%%

% r = rotate([0 0 -90])
r = [
    0.0000    1.0000         0         0
   -1.0000    0.0000         0         0
         0         0    1.0000         0
         0         0         0    1.0000
];

elec_ctf = ft_transform_geometry(r, elec); % rotate 90 degrees
elec_ctf.coordsys = 'ctf';

skin_ctf = ft_transform_geometry(r, skin); % rotate 90 degrees
skin_ctf.coordsys = 'ctf';

grad_ctf = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds'), 'senstype', 'meg');

%%


cfg = [];
cfg.template.headshape = skin_ctf;
cfg.template.headshapestyle = 'surface';

cfg.individual.elec = elec_ctf;
out = ft_interactiverealign(cfg);

%%

cfg = [];
cfg.template.headshape = skin_ctf;
cfg.template.headshapestyle = 'surface';

cfg.individual.grad = grad_ctf;
cfg.individual.gradstyle = {'coilshape', 'point'};
out = ft_interactiverealign(cfg);

cfg.individual.grad = grad_ctf;
cfg.individual.gradstyle = {'coilshape', 'circle'};
out = ft_interactiverealign(cfg);

