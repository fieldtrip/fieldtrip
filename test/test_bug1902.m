function test_bug1902

% WALLTIME 02:00:00
% MEM 3gb

% TEST ft_volumesegment ft_prepare_sourcemodel

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/nifti/single_subj_T1.nii'));
mri.coordsys = 'spm';

cfg=[];
tpm = ft_volumesegment(cfg,mri);  % gray, white, csf tissue prob. map

cfg=[];
cfg.mri = tpm;
grid01 = ft_prepare_sourcemodel(cfg,tpm);

cfg=[];
cfg.mri = tpm;
grid02 = ft_prepare_sourcemodel(cfg,tpm);

cfg=[];
cfg.mri = tpm;
grid03 = ft_prepare_sourcemodel(cfg,tpm);

grid01 = rmfield(grid01, 'cfg');
grid02 = rmfield(grid02, 'cfg');
grid03 = rmfield(grid03, 'cfg');

isequal(grid01,grid02);
isequal(grid01,grid03);
