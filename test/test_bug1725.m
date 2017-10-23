function test_bug1725

% MEM 2000mb
% WALLTIME 00:10:00

% TEST ft_read_atlas ft_prepare_atlas

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1725/afni'));

% try to read a collection of atlasses
% the first one is from AFNI
filename = 'TTatlas+tlrc.BRIK';
atlas = ft_read_atlas(filename);

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1725/wfu_pickatlas'));

% these ones are from http://fmri.wfubmc.edu/
filename = {
  % 'MNI_atlas_templates/MNI_T1.nii'
  'MNI_atlas_templates/TD_brodmann.nii'
  'MNI_atlas_templates/TD_hemisphere.nii'
  'MNI_atlas_templates/TD_label.nii'
  % 'MNI_atlas_templates/TD_label_extended_modified.nii'
  'MNI_atlas_templates/TD_lobe.nii'
  'MNI_atlas_templates/TD_type.nii'
  'MNI_atlas_templates/aal_MNI_V4.nii'
  'MNI_atlas_templates/atlas116.nii'
  'MNI_atlas_templates/atlas71.nii'
  % 'TD-ICBM_MNI_atlas_templates/MNI_T1.nii'
  'TD-ICBM_MNI_atlas_templates/TD_brodmann.nii'
  'TD-ICBM_MNI_atlas_templates/TD_hemisphere.nii'
  'TD-ICBM_MNI_atlas_templates/TD_label.nii'
  'TD-ICBM_MNI_atlas_templates/TD_lobe.nii'
  'TD-ICBM_MNI_atlas_templates/TD_type.nii'
  % 'mouse_atlas_templates/c57_15_M_Normal_age12.nii'
  'mouse_atlas_templates/mouseLabels.nii'
  'rhesus_atlas_templates/SubCortical.nii'
  % 'rhesus_atlas_templates/template.nii'
  % 'rhesus_atlas_templates/templateSkull.nii'
  'rhesus_atlas_templates/template_hardseg.nii'
  'rhesus_atlas_templates/template_parc.nii'
  % 'vervet_atlas_templates/Vervet_T1_Template_WholeHead.nii'
  'vervet_atlas_templates/vervetLabels.nii'
  'vervet_atlas_templates/vervetType.nii'
  };

for i=1:length(filename)%14:length(filename)
  % just try to read it
  disp(i);
  disp(filename{i});
  atlas = ft_read_atlas(filename{i});
end

% plot the template MRI together with the atlas
mri   = ft_read_mri('MNI_atlas_templates/MNI_T1.nii');
atlas = ft_read_atlas('MNI_atlas_templates/aal_MNI_V4.nii'); % FIXME this is how it should be

% FIXME is there a need to test the old to-be-deprecated functionality?
% atlas = ft_prepare_atlas('MNI_atlas_templates/aal_MNI_V4.nii'); % use the old representation

if false
  cfg = [];
  cfg.atlas = atlas;
  cfg.interactive = 'yes';
  ft_sourceplot(cfg, mri);
end

