function test_pull360

% MEM 2gb
% WALLTIME 00:10:00

mri				= ft_read_mri('single_subj_T1_1mm.nii', 'datatype', 'nifti');
mri.transform	= [-1,0,0,95.7;0,1,0,-122.5;0,0,1,-132.5;0,0,0,1]; % this will require intervention by getinside in ft_prepare_headmodels

cfg                         = [];
cfg.mri                     = mri;
cfg.sourcemodel.resolution         = 10;       % in mm
cfg.sourcemodel.unit               = 'mm';
ft_prepare_sourcemodel(cfg);
