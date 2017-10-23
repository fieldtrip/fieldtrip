function test_bug2986

% WALLTIME 00:20:00
% MEM 1500mb

% TEST ft_volumerealign ft_volumereslice

load standard_mri

% convert to 'cm' to show what happens
mri = ft_convert_units(mri, 'cm');

% reslice to standard (improves solutions) 
mri =  ft_volumereslice([], mri);

% load headshape
load(fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test'),'bug2986.mat'));

% for quick check feed in fiducial positions
fiducial = [];
fiducial.nas = [131 215 85];
fiducial.lpa =  [41 108 84];
fiducial.rpa = [213 107 80];
fiducial.zpoint = [131 125 195];

% coarse realignment based on fiducials first
cfg          = [];
cfg.method   = 'fiducial';
cfg.fiducial = fiducial;
cfg.coordsys = 'neuromag';
mri          = ft_volumerealign(cfg,mri);

% make automatic coregistration based on headshape
% the coregistration looks fine (depending on how well you defined the
% fiducials of course) - don't change anything in the interactive menu to
% avoid that that's the problem
cfg           = [];
cfg.method    = 'headshape';
cfg.headshape.headshape = shape;
cfg.headshape.interactive    = 'no';
cfg.headshape.icp = 'yes';
cfg.coordsys  = 'neuromag';
mri_aligned   = ft_volumerealign(cfg, mri);

cfg.headshape.shape = ft_convert_units(shape, 'm');
mri_aligned2  = ft_volumerealign(cfg, mri);
mri_aligned3  = ft_volumerealign(cfg, ft_convert_units(mri, 'mm'));


%% check based on fiducial + headshape alignment
% segment the scalp to check the alignment
cfg        = [];
cfg.output = 'scalp';
cfg.smooth = 2;
seg_align  = ft_volumesegment(cfg, mri_aligned);
seg_align2  = ft_volumesegment(cfg, mri_aligned2);
seg_align3  = ft_volumesegment(cfg, mri_aligned3);

% build headmodel
cfg             = [];
cfg.method      = 'singleshell';
cfg.numvertices = 2000;
hdm    = ft_prepare_headmodel(cfg, seg_align);
hdm2   = ft_prepare_headmodel(cfg, seg_align2);
hdm3   = ft_prepare_headmodel(cfg, seg_align3);

% now plot
figure; hold on
ft_plot_vol(ft_convert_units(hdm, 'cm'),'edgecolor','none','facecolor','w');
ft_plot_headshape(shape);

figure; hold on
ft_plot_vol(ft_convert_units(hdm2, 'cm'),'edgecolor','none','facecolor','w');
ft_plot_headshape(shape);

figure; hold on
ft_plot_vol(ft_convert_units(hdm3, 'cm'),'edgecolor','none','facecolor','w');
ft_plot_headshape(shape);

