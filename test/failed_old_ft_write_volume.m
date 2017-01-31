function failed_old_ft_write_volume

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_old_ft_write_volume

% this script tests whether ft_volumewrite and ft_write_volume works

% create some data
ft_hastoolbox('spm8',1);
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/external/spm8/templates/T1.nii'));

% create 'blob' in occipital cortex
tmp = zeros(mri.dim);
tmp(46,20,39) = 10;
spm_smooth(tmp, tmp, 10);
mri.pow = tmp;

% visualize it using fieldtrip
cfg = [];
cfg.funparameter = 'pow';
figure;ft_sourceplot(cfg, mri);

% write the functional volume to disk
fname = tempname;

% write to file
% read in the files again
cfg = [];
cfg.filename  = fname;
cfg.parameter = 'pow';

restoredefaultpath
cd ~jansch/matlab/fieldtrip
addpath(pwd);
ft_defaults
ft_hastoolbox('spm2', 1);
cfg.filetype  = 'analyze_spm';
ft_volumewrite(cfg, mri);
pow1 = ft_read_mri([fname,'.img']);

restoredefaultpath
cd ~jansch/matlab/fieldtrip
addpath(pwd);
ft_defaults
ft_hastoolbox('spm8', 1); % ensure using spm8
cfg.filetype  = 'nifti';
ft_volumewrite(cfg, mri);
pow2 = ft_read_mri([fname,'.nii']);

all(pow1.anatomy(:)==pow2.anatomy(:))
all(pow1.transform(:)==pow2.transform(:))

% both seem to be the same

mri.pow1 = pow1.anatomy;
mri.pow2 = pow2.anatomy;

cfg = [];
cfg.funparameter = 'pow1';
figure;ft_sourceplot(cfg, mri);
cfg.funparameter = 'pow2';
figure;ft_sourceplot(cfg, mri);

delete([fname,'.nii']);
delete([fname,'.img']);
delete([fname,'.hdr']);
