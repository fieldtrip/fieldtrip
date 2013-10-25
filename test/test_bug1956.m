function test_bug1902

% MEM 1gb
% WALLTIME 00:07:42

% TEST test_bug1902 ft_prepare_sourcemodel volumesmooth

mri = ft_read_mri('/home/common/matlab/fieldtrip/data/test/latest/mri/nifti/single_subj_T1.nii');

cfg=[];
cfg.coordsys = 'ctf';
tpm = ft_volumesegment(cfg,mri);  % gray, white, csf tissue prob. map
cfg=[];
cfg.mri=tpm;
grid01 = ft_prepare_sourcemodel(cfg,tpm);

cfg=[];
cfg.mri=tpm;
grid02 = ft_prepare_sourcemodel(cfg,tpm);

cfg=[];
cfg.mri=tpm;
grid03 = ft_prepare_sourcemodel(cfg,tpm);

grid01 = rmfield(grid01,'cfg');
grid02 = rmfield(grid02, 'cfg');

isequal(grid01,grid02)


