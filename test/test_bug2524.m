function test_bug2524

% WALLTIME 00:10:00
% MEM 1024mb

% TEST ft_write_mri

data = zeros(5,5,5);
trans = diag([5 5 5 1]);

tmp = tempname;
ft_write_mri([tmp,'.nii'],data,'transform',trans,'dataformat','nifti');
mri = ft_read_mri([tmp,'.nii']);

assert(mri.transform(1,1)==5);

% use a command line tool from FSL
[r,s] = system('which fslhd');
if r==0
  disp(r)
  disp(s)
else
  % use hard-coded FSL version
  s = '/opt/fsl/5.0.9/bin/fslhd';
end
[r,s] = system([deblank(s) ' ' tmp '.nii | grep pixdim']);
disp(r)
disp(s)

assert(numel(strfind(s,'5.0'))==3);

delete([tmp,'.nii']);
