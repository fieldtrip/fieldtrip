function test_bug2909

% WALLTIME 00:10:00
% MEM 3gb

% TEST test_bug2865
% TEST ft_read_cifti

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2909'));

%%
% this is the data mentioned by JM on http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2865
filename1 = 'RSN-networks.32k_fs_LR.dlabel.nii';
filename2 = 'RSN-networks.32k_fs_LR.dlabel.nii';
surfL = 'Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii';		
surfR = 'Q1-Q6_R440.R.midthickness.32k_fs_LR.surf.gii';

cii1 = ft_read_cifti(filename1, 'cortexleft', surfL, 'cortexright', surfR, 'mapname', 'field');
cii2 = ft_read_cifti(filename2, 'cortexleft', surfL, 'cortexright', surfR, 'mapname', 'array');

% this is the file I got from Mike Harms over email
filename3 = 'rfMRI_REST1_LR_Atlas_stats_ALLCOL_MEAN.dscalar.nii';
filename4 = 'rfMRI_REST1_LR_Atlas_stats_ALLCOL_MEAN.dscalar.nii';

cii3 = ft_read_cifti(filename3, 'debug', true, 'maplabel', 'field');
cii4 = ft_read_cifti(filename4, 'debug', true, 'maplabel', 'array');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try writing them and subsequently reading them again

filename1b = tempname;
ft_write_cifti(filename1b, cii1, 'parameter', 'rsn_consensus_communities___pcn11__power_neuron11');
d = dir([filename1b '*.nii']);
filename1b = fullfile(fileparts(filename1b), d.name); % determine the name of the cii file
cii1b = ft_read_cifti(filename1b, 'cortexleft', surfL, 'cortexright', surfR);
assert(isequal(cii1.rsn_consensus_communities___pcn11__power_neuron11, cii1b.dscalar));

%% this one does not work yet
% filename2b = tempname;
% ft_write_cifti(filename2b, cii2, 'parameter', 'x32k_fs_lr');
% d = dir([filename2b '*.nii']);
% filename1b = fullfile(fileparts(filename2b), d.name); % determine the name of the cii file
% cii2b = ft_read_cifti(filename1b, 'cortexleft', surfL, 'cortexright', surfR);






