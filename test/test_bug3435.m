function test_bug3435

% WALLTIME 00:20:00
% MEM 3gb

% TEST prepare_mesh_tetrahedral.m


mri = ft_read_mri('/home/common/matlab/fieldtrip/template/anatomy/single_subj_T1.nii');

cfg = [];
mri.coordsys = 'ras';
seg_mri = ft_volumesegment(cfg,mri);

mesh = prepare_mesh_tetrahedral([],seg_mri);
