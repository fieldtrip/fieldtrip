function test_bug2741

% WALLTIME 00:10:00
% MEM 1500mb

% TEST ft_read_cifti ft_write_cifti

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2741'));

cii1 = ft_read_cifti('zstat1.dtseries.nii', 'cortexleft', '100307.L.midthickness.32k_fs_LR.surf.gii', 'cortexright', '100307.R.midthickness.32k_fs_LR.surf.gii');

ft_write_cifti('zstat2', cii1, 'parameter', 'dtseries');

cii2 = ft_read_cifti('zstat2.dtseries.nii'); % surfaces will be read automatically

cii1 = rmfield(cii1, 'hdr');
cii2 = rmfield(cii2, 'hdr');

assert(isequaln(cii1, cii2));
