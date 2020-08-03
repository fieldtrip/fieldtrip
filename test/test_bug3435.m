function test_bug3435

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY ft_prepare_mesh prepare_mesh_tetrahedral

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/template/anatomy/single_subj_T1.nii'));
mri.coordsys = 'mni';

cfg = [];
cfg.downsample = 2;
cfg.output = {'scalp', 'skull', 'csf', 'gray', 'white'};
mri_segmented = ft_volumesegment(cfg, mri);

cfg = [];
cfg.method = 'tetrahedral';
mesh = ft_prepare_mesh(cfg, mri_segmented);

figure
ft_plot_mesh(mesh, 'surfaceonly', 1, 'facecolor', 'w');
