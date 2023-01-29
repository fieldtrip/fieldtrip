function test_ft_volumewrite

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_volumewrite

[ftver, ftpath] = ft_version;
mrifilename = fullfile(ftpath, 'template', 'headmodel', 'standard_mri.mat');

mri = ft_read_mri(mrifilename);

cfg = [];
cfg.parameter = 'anatomy';
cfg.filename = tempname;
cfg.filetype = 'nifti_spm';
ft_volumewrite(cfg, mri)
delete([cfg.filename '*'])
