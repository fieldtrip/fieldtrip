function test_bug1832
% TEST test_bug1832
% TEST ft_read_mri ft_volumesegment ft_prepare_headmodel ft_prepare_sourcemodel
%
% for the warp template example script:
% 1) make test script that checks on the units in these objects (with assert)
% 2) also go through the rest of the example script and check on units


% NOTE: the path to the template file is user-specific
if isunix
  template = ft_read_mri('home/common/matlab/fieldtrip/external/spm8/templates/T1.nii');
elseif ispc
  template = ft_read_mri('H:/common/matlab/fieldtrip/external/spm8/templates/T1.nii');
end
template.coordsys = 'spm'; % so that FieldTrip knows how to interpret the coordinate system

% segment the template brain and construct a volume conduction model (i.e. head model): this is needed
% for the inside/outside detection of voxels.
cfg          = [];
template_seg = ft_volumesegment(cfg, template);
assert(isfield(template_seg, 'unit'), 'unit field in seg missing')

cfg          = [];
cfg.method   = 'singleshell';
template_vol = ft_prepare_headmodel(cfg, template_seg);
assert(isfield(template_vol, 'unit'), 'unit field in vol missing')
% construct the dipole grid in the template brain coordinates
% the source units are in cm
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg = [];
cfg.grid.xgrid  = -20:1:20;
cfg.grid.ygrid  = -20:1:20;
cfg.grid.zgrid  = -20:1:20;
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.vol        = template_vol;
template_grid  = ft_prepare_sourcemodel(cfg);
assert(isfield(template_grid, 'unit'), 'unit field in grid missing')

% check more here