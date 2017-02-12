function test_bug2032

% WALLTIME 00:20:00
% MEM 2gb

% read each of the atlasses and ensure that it specifies the units and coordsys

[ftver, ftpath] = ft_version;

%%

filename = {
  'template/atlas/aal/ROI_MNI_V4.nii'
  'template/atlas/afni/TTatlas+tlrc.BRIK'
  'template/atlas/brainweb/brainweb_discrete.mat'
  'template/atlas/brainweb/brainweb_fuzzy.mat'
  'template/atlas/spm_anatomy/AllAreas_v17.img'
  'template/atlas/spm_anatomy/AllAreas_v18.img'
  'template/atlas/vtpm/vtpm.mat'
  };

for i=1:numel(filename)
  atlas = ft_read_atlas(fullfile(ftpath, filename{i}));
  assert(isfield(atlas, 'unit'), sprintf('unit missing in %s', filename{i}));
  assert(isfield(atlas, 'coordsys'), sprintf('coordsys missing in %s', filename{i}));
end

%%

% the inflated ones do not have a coordsys
%  'template/anatomy/surface_inflated_both.mat'
%  'template/anatomy/surface_inflated_both_caret.mat'
%  'template/anatomy/surface_inflated_left.mat'
%  'template/anatomy/surface_inflated_left_caret.mat'
%  'template/anatomy/surface_inflated_right.mat'
%  'template/anatomy/surface_inflated_right_caret.mat'


filename = {
  'template/anatomy/surface_pial_both.mat'
  'template/anatomy/surface_pial_left.mat'
  'template/anatomy/surface_pial_right.mat'
  'template/anatomy/surface_white_both.mat'
  'template/anatomy/surface_white_left.mat'
  'template/anatomy/surface_white_right.mat'
  };

for i=1:numel(filename)
  surface = ft_read_headshape(fullfile(ftpath, filename{i}));
  assert(isfield(surface, 'unit'), sprintf('unit missing in %s', filename{i}));
  assert(isfield(surface, 'coordsys'), sprintf('coordsys missing in %s', filename{i}));
end
