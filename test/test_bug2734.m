function test_bug2734

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_cifti


cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2734'));

cii = ft_read_cifti('tstat1.dtseries.nii');

assert(size(cii.pos,1)==96854);
assert(numel(cii.dtseries)==96854);
assert(numel(unique(cii.brainstructurelabel))==max(cii.brainstructure));
