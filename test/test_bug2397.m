function test_bug2397

% MEM 2000mb
% WALLTIME 00:10:00

% TEST ft_prepare_mesh

% this pertains to the OpenMEEG pipeline, see http://bugzilla.fcdonders.nl/show_bug.cgi?id=2397

% use the same test data as test_bug2396
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2396/segmentation/brainsuite/nobias_KR.nii'));
seg = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2396/segmentation/brainsuite/nobias_KR.label.nii'));

% dress up the segmentation
seg.tissue = seg.anatomy;
seg = rmfield(seg, 'anatomy');
seg.tissuelabel = {};
seg.tissuelabel{19} = 'brain';
seg.tissuelabel{18} = 'csf';
seg.tissuelabel{17} = 'skull';
seg.tissuelabel{16} = 'scalp';

cfg = [];
cfg.method = 'iso2mesh';
cfg.numvertices = 500; % use the same value for each mesh
mesh = ft_prepare_mesh(cfg, seg); % it is ugly, but runs through

cfg = [];
cfg.method = 'iso2mesh';
cfg.tissue = {'scalp'  'skull'  'csf'  'brain'};
cfg.numvertices = [1 2 3 3]*1000; % the inside meshes are more refined
mesh = ft_prepare_mesh(cfg, seg); % again ugly

cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = {'scalp'  'skull'  'csf'  'brain'};
cfg.numvertices = [1 2 3 3]*1000; % the inside meshes are more refined
mesh = ft_prepare_mesh(cfg, seg); % much better, although not correct
