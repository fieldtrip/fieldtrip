function test_issue1459

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_checkdata ft_datatype_segmentation ft_read_atlas

[ftver, ftpath] = ft_version;

%%

filename = {
  fullfile(ftpath, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii')
  fullfile(ftpath, 'template', 'atlas', 'afni', 'TTatlas+tlrc.HEAD')
  fullfile(ftpath, 'template', 'atlas', 'brainweb', 'brainweb_discrete.mat')
  fullfile(ftpath, 'template', 'atlas', 'brainweb', 'brainweb_fuzzy.mat')
  fullfile(ftpath, 'template', 'atlas', 'brainnetome', 'BNA_MPM_thr25_1.25mm.nii')
  fullfile(ftpath, 'template', 'atlas', 'spm_anatomy', 'AllAreas_v17.img')
  fullfile(ftpath, 'template', 'atlas', 'spm_anatomy', 'AllAreas_v18.img')
  fullfile(ftpath, 'template', 'atlas', 'vtpm', 'vtpm.mat')
  fullfile(ftpath, 'template', 'atlas', 'yeo', 'Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii')
  fullfile(ftpath, 'template', 'atlas', 'yeo', 'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii')
  };

%%

for i=1:numel(filename)
  disp('-------------------------------------------------------------------------------------');
  disp(filename{i});
  disp('-------------------------------------------------------------------------------------');
  atlas = ft_read_atlas(filename{i});
  ft_checkdata(atlas, 'feedback', 'yes');
end
