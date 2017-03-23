function test_pull360
% TEST ft_prepare_sourcemodel

% WALLTIME 00:10:00
% MEM 500mb

mri				= ft_read_mri('single_subj_T1_1mm.nii', 'datatype', 'nifti');
mri.transform	= [-1,0,0,95.7;0,1,0,-122.5;0,0,1,-132.5;0,0,0,1]; % this will require intervention by getinside in ft_prepare_headmodels

cfg                         = [];
cfg.mri                     = mri;
cfg.grid.resolution         = 10;       % in mm
cfg.grid.unit               = 'mm';
ft_prepare_sourcemodel(cfg);
