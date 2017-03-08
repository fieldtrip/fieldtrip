function test_bug1828

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_convert_coordsys
% TEST align_ctf2spm

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1828'));

acvox = [89 135 125]; % voxel that is identified as ~ AC, i.e. the origin of the coordinate system

% make sure a version of SPM is on the path
ft_hastoolbox('SPM8',1);

mri0 = ft_convert_coordsys(mri, 'spm', 0);
mri1 = ft_convert_coordsys(mri, 'spm', 1);
mri2 = ft_convert_coordsys(mri, 'spm', 2);

acvox0 = ft_warp_apply(inv(mri0.transform), [0 0 0]);
acvox1 = ft_warp_apply(inv(mri1.transform), [0 0 0]);
acvox2 = ft_warp_apply(inv(mri2.transform), [0 0 0]);

if 0,
  figure;hold on;
  imagesc(mri.anatomy(:,:,round(acvox(3))));colormap gray;
  plot(round(acvox(2)),round(acvox(1)),'ro');
  figure;hold on;
  imagesc(mri.anatomy(:,:,round(acvox0(3))));colormap gray;
  plot(round(acvox0(2)),round(acvox0(1)),'yo');
  figure;hold on;
  imagesc(mri.anatomy(:,:,round(acvox1(3))));colormap gray;
  plot(round(acvox1(2)),round(acvox1(1)),'yo');
  figure;hold on;
  imagesc(mri.anatomy(:,:,round(acvox2(3))));colormap gray;
  plot(round(acvox2(2)),round(acvox2(1)),'yo');
end
  
err0 = sqrt(sum((acvox-acvox0).^2));
err1 = sqrt(sum((acvox-acvox1).^2));
err2 = sqrt(sum((acvox-acvox2).^2));

assert(err2<err1);
assert(err2<err0);

