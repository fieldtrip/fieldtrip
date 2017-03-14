function test_bug1832

% MEM 1500mb
% WALLTIME 00:20:00

% TEST ft_read_mri ft_volumesegment ft_prepare_headmodel ft_prepare_sourcemodel

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% for the warp template example script:
% 1) make test script that checks on the units in these objects (with assert)
% 2) also go through the rest of the example script and check on units

template = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/external/spm8/templates/T1.nii'));
template.coordsys = 'spm'; % so that FieldTrip knows how to interpret the coordinate system

% segment the template brain and construct a volume conduction model (i.e. head model): this is needed
% for the inside/outside detection of voxels.
cfg          = [];
template_seg = ft_volumesegment(cfg, template);
assert(isfield(template_seg, 'unit'), 'unit field in seg missing')

cfg          = [];
cfg.method   = 'singlesphere';
template_vol = ft_prepare_headmodel(cfg, template_seg);
assert(isfield(template_vol, 'unit'), 'unit field in vol missing')
% construct the dipole grid in the template brain coordinates
% the source units are in mm
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg = [];
cfg.grid.xgrid  = -200:10:200;
cfg.grid.ygrid  = -200:10:200;
cfg.grid.zgrid  = -200:10:200;
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.vol        = template_vol;
template_grid  = ft_prepare_sourcemodel(cfg);
assert(isfield(template_grid, 'unit'), 'unit field in grid missing')

% prepare leadfield
cfg         = [];
cfg.grid    = ft_convert_units(template_grid, 'cm');
cfg.vol     = template_vol;
cfg.channel = 'EEG';
cfg.elec = ft_read_sens('standard_1020.elc');
template_grid_lf     = ft_prepare_leadfield(cfg);
assert(isfield(template_grid_lf, 'unit'), 'unit field in leadfield missing')

template_grid_lf = ft_convert_units(template_grid_lf, template_grid.unit);
assert(all(template_grid_lf.pos(:) == template_grid.pos(:)), 'leadfield positions differ from sourcemodel position');
