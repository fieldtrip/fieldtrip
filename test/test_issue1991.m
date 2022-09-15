function test_issue1991

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_read_mri ft_hastoolbox ft_filetype

[ft_ver, ft_path] = ft_version;

file0 = fullfile(ft_path, 'template', 'anatomy', 'single_subj_T1_1mm.nii');
file1 = dccnpath('/home/common/matlab/fieldtrip/data/test/original/jnifti/samples/colin27/colin27_zlib.jnii');
file2 = dccnpath('/home/common/matlab/fieldtrip/data/test/original/jnifti/samples/colin27/colin27_zlib.bnii');

%%

mri0 = ft_read_mri(file0);
mri1 = ft_read_mri(file1);
mri2 = ft_read_mri(file2);

%%

% do not compare the anatomy in the samples to the original, as that is obviously different
% it appears that the jnifti/samples represent a segmentation rather than the original MRI
fn = {'dim', 'unit', 'transform'};

for i=1:length(fn)
  assert(isequal(mri0.(fn{i}), mri1.(fn{i})), ['difference in jnii field ' fn{i}]);
  assert(isequal(mri0.(fn{i}), mri2.(fn{i})), ['difference in bnii field ' fn{i}]);
end

%%

assert(isequal(mri1.anatomy, mri2.anatomy), 'difference in anatomy')

%%

cfg = [];
figure; ft_sourceplot(cfg, mri0);
figure; ft_sourceplot(cfg, mri1);
figure; ft_sourceplot(cfg, mri2);
