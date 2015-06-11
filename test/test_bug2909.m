function test_bug2865

% WALLTIME 0:10:00
% MEM 3gb

% TEST test_bug2865
% TEST ft_read_cifti

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2909'));

% this is the data mentioned by JM on http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2865
filename = 'RSN-networks.32k_fs_LR.dlabel.nii';
surfL = 'Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii';		
surfR = 'Q1-Q6_R440.R.midthickness.32k_fs_LR.surf.gii';

cii = ft_read_cifti(filename, 'cortexleft', surfL, 'cortexright', surfR);

% this is the file I got from Mike Harms over email
filename = 'rfMRI_REST1_LR_Atlas_stats_ALLCOL_MEAN.dscalar.nii';

cii1 = ft_read_cifti(filename, 'debug', true, 'maplabel', 'field');
cii2 = ft_read_cifti(filename, 'debug', true, 'maplabel', 'array');

