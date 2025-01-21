function test_ft_determine_coordsys

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY
% DATA private

mrifile = dccnpath('/project/3031000.02/test/latest/mri/nifti/single_subj_T1.nii');
mri     = ft_read_mri(mrifile);

ft_determine_coordsys(mri,'interactive','no')
close

ft_determine_coordsys(mri,'clim',[0 .1],'interactive','no')
close

