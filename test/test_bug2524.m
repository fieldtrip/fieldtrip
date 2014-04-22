function test_bug2524

% WALLTIME 00:10:00
% MEM 1024mb

% TEST test_bug2524
% TEST ft_write_mri

data = zeros(5,5,5);
trans = diag([5 5 5 1]);

tmp = tempname;
ft_write_mri([tmp,'.nii'],data,'transform',trans,'dataformat','nifti');
mri = ft_read_mri([tmp,'.nii']);

assert(mri.transform(1,1)==5);
[r,s] = system(['fslhd ',tmp,'.nii | grep pixdim']);
assert(numel(strfind(s,'5.0'))==3);

delete([tmp,'.nii']);
