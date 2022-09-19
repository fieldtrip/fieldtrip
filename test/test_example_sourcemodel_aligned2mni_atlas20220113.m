function test_example_sourcemodel_aligned2mni_atlas

% MEM 4gb
% WALLTIME 00:10:00

%
%% Create brain atlas based MNI-aligned grids in individual head-space
%
% When combining the source-level data of multiple subjects this data is typically first interpolated (e.g., using **[ft_sourceinterpolate](https://github.com/fieldtrip/fieldtrip/blob/release/ft_sourceinterpolate.m)**) and then spatially normalized to a template brain (e.g., using **[ft_volumenormalise](https://github.com/fieldtrip/fieldtrip/blob/release/ft_volumenormalise.m)**). However it is also possible to define the source-reconstruction grid for each individual subject in such a way that all these grids are already aligned in MNI-space. The combination or statistic of source-level data across subjects can then directly be computed within the source-structure without the need to interpolate and normalize each volume. In addition, the position of the grid points can be chosen in a way that only locations corresponding to particular brain areas (parcels) are included.
%
% The procedure is as follows. First, a template grid is computed on the basis of the standard head model located in the FieldTrip template directory. Subsequently, a brain atlas is loaded using **[ft_read_atlas](https://github.com/fieldtrip/fieldtrip/blob/release/fileio/ft_read_atlas.m)**. On the basis on the brain atlas and the template grid it is possible to create a binary mask of all locations in the template grid that match atlas locations using **[ft_volumelookup](https://github.com/fieldtrip/fieldtrip/blob/release/ft_volumelookup.m)**. Finally the mask is used to determine which location shall be defined as 'inside' the head model.
%
[ftver, ftpath] = ft_version;
load(fullfile(ftpath, 'template/headmodel/standard_singleshell.mat'));

cfg = [];
cfg.xgrid  = -20:1:20;
cfg.ygrid  = -20:1:20;
cfg.zgrid  = -20:1:20;
cfg.unit   = 'cm';
cfg.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.headmodel   = vol;
template_grid   = ft_prepare_sourcemodel(cfg);
template_grid.coordsys = 'mni';

template_grid = ft_convert_units(template_grid,'cm');

figure;
hold on
ft_plot_mesh(template_grid.pos(template_grid.inside,:));
ft_plot_headmodel(vol,  'facecolor', 'cortex', 'edgecolor', 'none');
ft_plot_axes(vol);
alpha 0.5
camlight

% Read the atlas, convert to units of cm and create the binary mask.
%
atlas = ft_read_atlas(fullfile(ftpath, 'template/atlas/aal/ROI_MNI_V4.nii'));

atlas = ft_convert_units(atlas,'cm');

cfg = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
mask           = ft_volumelookup(cfg, template_grid);

% Now we determine all indices of the binary mask to be considered as inside the head model. And plot the result. Note the missing dipole locations for example in the vicinity of the ventricles.
%
template_grid.inside = false(template_grid.dim);
template_grid.inside(mask==1) = true;

figure;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%
% Load the subject-specific MRI from [here](ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/salzburg/mri.mat) and inverse-warp the subject specific grid to the template grid.
%
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/salzburg/mri.mat'));

cfg                = [];
cfg.warpmni   = 'yes';
cfg.template  = template_grid;
cfg.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri;
sourcemodel        = ft_prepare_sourcemodel(cfg);

% Finally, you can load the subject-specific headmodel from [here](ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/salzburg/hdm.mat) and check the result with the following code.
%
close all

hdm = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/salzburg/hdm.mat'));

hdm         = ft_convert_units(hdm, 'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');

figure
hold on
ft_plot_headmodel(hdm,  'facecolor', 'cortex', 'edgecolor', 'none');
ft_plot_axes(hdm);
alpha 0.4  % make the surface transparent
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
view ([0 -90 0])
