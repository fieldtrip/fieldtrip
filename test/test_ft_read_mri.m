function test_ft_read_mri

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_read_mri

%%

files = {
  'afni/anat+orig.BRIK'
  'afni/anat+orig.HEAD'
  'analyze/oostenveld_r.hdr'
  'analyze/oostenveld_r.img'
  'ctf_v22/anonymous.mri'
  'ctf_v41/anonymous.mri'
  'dicom/19112010_JHORSCHIG.MR.FCDC_SEQUENCES_STANDARD_SEQUENCES.0002.0064.2010.11.19.12.08.01.265625.73005239.IMA'
  'freesurfer/T1.mgz'
  'minc/single_subj_T1.mnc'
  'nifti/single_subj_T1.nii'
  'neuromag/scans/mri.fif'
  };
  %'neuromag/slices/MR1.3.12.2.1107.5.2.32.35204.2008010817494647729256323'

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri');

for k = 1:numel(files)
  filename = dccnpath(fullfile(datadir,files{k}));
  ft_read_mri(filename);
end

%%

[v, p] = ft_version;
filename = fullfile(p, 'template', 'anatomy', 'single_subj_T1.nii');

mri0 = ft_read_mri(filename, 'dataformat', 'nifti');

rmpath(fileparts(which('spm')))
mri1 = ft_read_mri(filename, 'dataformat', 'nifti_spm', 'spmversion', 'spm8');

rmpath(fileparts(which('spm')))
mri2 = ft_read_mri(filename, 'dataformat', 'nifti_spm', 'spmversion', 'spm12');

% they should all be the same
assert(isequal(mri0.anatomy, mri1.anatomy));
assert(isequal(mri0.anatomy, mri2.anatomy));

% they should all be the same
assert(isequal(mri0.transform, mri1.transform));
assert(isequal(mri0.transform, mri2.transform));
