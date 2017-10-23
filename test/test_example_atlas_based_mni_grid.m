function test_example_atlas_based_mni_grid

% WALLTIME 00:20:00
% MEM 3gb

% http://www.fieldtriptoolbox.org/example/create_single-subject_grids_in_individual_head_space_that_are_all_aligned_in_brain_atlas_based_mni_space

%%

[ftver, ftpath] = ft_version;

load(fullfile(ftpath, 'template/headmodel/standard_singleshell.mat'));

%%

cfg = [];
cfg.grid.xgrid  = -20:1:20;
cfg.grid.ygrid  = -20:1:20;
cfg.grid.zgrid  = -20:1:20;
cfg.grid.unit   = 'cm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.headmodel   = vol;
template_grid   = ft_prepare_sourcemodel(cfg);

template_grid = ft_convert_units(template_grid,'cm');

figure;
hold on
ft_plot_mesh(template_grid.pos(template_grid.inside,:));
ft_plot_vol(vol,  'facecolor', 'cortex', 'edgecolor', 'none');
ft_plot_axes(vol);
alpha 0.5
camlight

%%

atlas = ft_read_atlas(fullfile(ftpath, 'template/atlas/aal/ROI_MNI_V4.nii'));

atlas = ft_convert_units(atlas, 'cm');

cfg = [];
cfg.atlas       = atlas;
cfg.roi         = atlas.tissuelabel; % here you can also specify a single label, i.e. single ROI
cfg.inputcoord  = 'mni';
mask            = ft_volumelookup(cfg,template_grid);

%%

template_grid.inside = false(template_grid.dim);
template_grid.inside(mask==1) = true;

figure;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%%

% data from ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/salzburg/mri.mat

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/salzburg/mri.mat'));

cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri;
sourcemodel        = ft_prepare_sourcemodel(cfg);

%%

close all

hdm = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/salzburg/hdm.mat'));

hdm         = ft_convert_units(hdm, 'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');

figure
hold on
ft_plot_vol(hdm,  'facecolor', 'cortex', 'edgecolor', 'none');
ft_plot_axes(hdm);
alpha 0.4  % make the surface transparent
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
view ([0 -90 0])

