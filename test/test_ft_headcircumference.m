function test_ft_headcircumference

% WALLTIME 00:10:00
% DEPENDENCY ft_headcircumference ft_volumesegment ft_prepare_mesh
% DATA no

filename = dccnpath('/project/3031000.02/external/download/test/ctf/Subject01.mri');
mri = ft_read_mri(filename);

%%

cfg = [];
cfg.spmversion = 'spm12';
cfg.output = 'scalp';
mri_segmented = ft_volumesegment(cfg, mri);

%%

cfg = [];
cfg.tissue = 'scalp';
cfg.method = 'projectmesh';
mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);

%%

figure
ft_plot_mesh(mesh_scalp)

%%
% not really needed, but the code is largely shared

fid = [
  117.9935   -2.5456   -1.5713
   11.1552  -74.6459   -8.1709
   15.4973   76.9586   -1.3693
  -78.6502    2.5375   31.1344
  ];

cfg = [];
cfg.method = '1020';
cfg.fiducial.nas   = fid(1,:);
cfg.fiducial.ini   = fid(4,:);
cfg.fiducial.lpa   = fid(3,:);
cfg.fiducial.rpa   = fid(2,:);
elec1020 = ft_electrodeplacement(cfg, mesh_scalp);

%%
% test on the triangulated mesh

cfg = [];
cfg.fiducial.nas   = fid(1,:);
cfg.fiducial.ini   = fid(4,:);
cfg.fiducial.lpa   = fid(3,:);
cfg.fiducial.rpa   = fid(2,:);
circumference = ft_headcircumference(cfg, mesh_scalp);

% the MRI, segmentation and mesh are in mm
assert(circumference>500);
assert(circumference<700);

%%
% construct a headshape like the one from a Polhemus
sel = mesh_scalp.pos(:,3)>-20;
headshape = [];
headshape.pos = mesh_scalp.pos(sel, :);

figure
ft_plot_mesh(headshape.pos)

%%
% test on the headshape as structure

cfg = [];
cfg.fiducial.nas   = fid(1,:);
cfg.fiducial.ini   = fid(4,:);
cfg.fiducial.lpa   = fid(3,:);
cfg.fiducial.rpa   = fid(2,:);
circumference = ft_headcircumference(cfg, headshape);

assert(circumference>500);
assert(circumference<700);

%%
% test on the headshape as Nx3 array

cfg = [];
cfg.fiducial.nas   = fid(1,:);
cfg.fiducial.ini   = fid(4,:);
cfg.fiducial.lpa   = fid(3,:);
cfg.fiducial.rpa   = fid(2,:);
circumference = ft_headcircumference(cfg, headshape.pos);

assert(circumference>500);
assert(circumference<700);