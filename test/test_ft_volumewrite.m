function test_ft_volumewrite

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_volumewrite

mrifilename   = dccnpath('/home/common/matlab/fieldtrip/template/headmodel/standard_mri.mat');

mri = ft_read_mri(mrifilename);
cfg = [];
cfg.parameter = 'anatomy';
cfg.filename = 'tempmri';
cfg.filetype = 'nifti';
ft_volumewrite(cfg, mri)
delete tempmri.nii

