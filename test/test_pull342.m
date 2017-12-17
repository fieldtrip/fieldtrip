function test_pull342

% WALLTIME 00:10:00
% MEM 2gb

%% Setup

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/pull342'));

mri = ft_read_mri('IR29_Bext.nii.gz');
% elec = read_bioimage_mgrid('IR29_grid.mgrid');
elec = ft_read_sens('IR29_grid.mgrid');

%% Visualize

% plot it at the middle of the volume
location = ft_warp_apply(mri.transform, mri.dim/2);

figure
ft_plot_ortho(mri.anatomy, 'transform', mri.transform, 'location', location, 'style', 'intersect');
ft_plot_axes(mri);
ft_plot_sens(elec);

% The file reading works, but the figure is not what I would expect. The electrodes
% appear displaced in the positive x-direction, the MRI is displaced in the negative
% x-direction. Neither is correct with respect to the origin at [0,0,0].

