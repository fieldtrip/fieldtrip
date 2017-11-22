function inspect_bug1830

% WALLTIME 00:10:00
% MEM 1gb

%%

[v, p] = ft_version;

mri = ft_read_mri(fullfile(p, 'template', 'anatomy', 'single_subj_T1.nii'));
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

% move the electrodes a little bit away from the original
elec = ft_transform_geometry(t * r, elec);

% move the MRI a little bit away from the original
mri = ft_transform_geometry(t * r, mri);


% figure
% ft_plot_mesh(skin, 'facecolor', 'skin');
% ft_plot_sens(elec);
% alpha 0.5


%%

cfg = [];
cfg.template.headshape = skin;
cfg.template.headshapestyle = 'surface';

cfg.template.axes = 'yes';
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
cfg.method = 'interactive';
cfg.elec = ft_convert_units(elec_ctf, 'cm'); % can also be passed as 2nd input data element
cfg.headshape = skin_ctf;
elec_realigned = ft_electroderealign(cfg);

assert(isequal(elec_realigned.unit, 'cm')); % should not have changed


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

%%
% following the redesign and cleanup of ft_interactiverealign, the goal is to call
% interactiverealign from the other data-type specific realign functions:
%  - ft_electroderealign.m	(method=interactive)
%  - ft_volumerealign.m     (for the optional interactive alignment in method=headshape)
%  - ft_sensorrealign.m	    (has been deprecated and will not get any further updates)
%  - ft_surfacerealign.m    (should be renamed in ft_meshrealign)

cfg = [];
cfg.method = 'interactive';
cfg.elec = ft_convert_units(elec_ctf, 'cm'); % can also be passed as 2nd input data element
cfg.target = elec;
elec_realigned = ft_electroderealign(cfg);

assert(isequal(elec_realigned.unit, 'cm')); % should not have changed

%%

cfg = [];
cfg.method = 'headshape';
cfg.headshape.interactive = 'yes';
cfg.headshape.icp = 'yes';
cfg.headshape.headshape = skin;
cfg.coordsys = 'mni';
cfg.spmversion = 'spm12';
mri_realigned = ft_volumerealign(cfg, mri);

% only interactive final check
cfg.headshape.icp = 'no';
mri_realigned2 = ft_volumerealign(cfg, mri_realigned);

%%

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
out = ft_meshrealign(cfg, skin);

ft_determine_coordsys(skin, 'interactive', 'no')
ft_determine_coordsys(out, 'interactive', 'no')

%%

cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
cfg.fiducial.nas = [2.5504 84.2050 -43.6891];
cfg.fiducial.lpa = [-85.1794 -24.1196 -42.7575];
cfg.fiducial.rpa = [85.6140 -15.1561 -49.6673];
%     zpoint = [-8.4323 -19.5475 100.8559]
out = ft_meshrealign(cfg, skin);


ft_determine_coordsys(skin, 'interactive', 'no')
ft_determine_coordsys(out, 'interactive', 'no')


%%

cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
out = ft_meshrealign(cfg, skin);

ft_determine_coordsys(skin, 'interactive', 'no')
ft_determine_coordsys(out, 'interactive', 'no')


