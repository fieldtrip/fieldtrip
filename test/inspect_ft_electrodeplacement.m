function inspect_ft_electrodeplacement

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY

%%
% do some prepatations

filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.mri');
mri = ft_read_mri(filename);

cfg             = [];
cfg.output      = 'scalp' ;
segmentation    = ft_volumesegment(cfg, mri);

cfg             = [];
cfg.method      = 'isosurface';
cfg.tissue      = 'scalp' ;
headshape       = ft_prepare_mesh(cfg, segmentation);

label = {
  'LPA'
  'RPA'
  'NAS'
  'INI'
  'Nz'
  'Fp1'
  'Fpz'
  'Fp2'
  'AF9'
  'AF7'
  'AF5'
  'AF3'
  'AF1'
  'AFz'
  'AF2'
  'AF4'
  'AF6'
  'AF8'
  'AF10'
  'F9'
  'F7'
  'F5'
  'F3'
  'F1'
  'Fz'
  'F2'
  'F4'
  'F6'
  'F8'
  'F10'
  'FT9'
  'FT7'
  'FC5'
  'FC3'
  'FC1'
  'FCz'
  'FC2'
  'FC4'
  'FC6'
  'FT8'
  'FT10'
  'T9'
  'T7'
  'C5'
  'C3'
  'C1'
  'Cz'
  'C2'
  'C4'
  'C6'
  'T8'
  'T10'
  'TP9'
  'TP7'
  'CP5'
  'CP3'
  'CP1'
  'CPz'
  'CP2'
  'CP4'
  'CP6'
  'TP8'
  'TP10'
  'P9'
  'P7'
  'P5'
  'P3'
  'P1'
  'Pz'
  'P2'
  'P4'
  'P6'
  'P8'
  'P10'
  'PO9'
  'PO7'
  'PO5'
  'PO3'
  'PO1'
  'POz'
  'PO2'
  'PO4'
  'PO6'
  'PO8'
  'PO10'
  'O1'
  'Oz'
  'O2'
  };

%%
% manual clicking of electrodes on on anatomical MRI

cfg = [];
cfg.method = 'volume';
cfg.channel = label;
elec1 = ft_electrodeplacement(cfg, mri);

%%
% manual clicking of electrodes on scalp surface

cfg = [];
cfg.method = 'headshape';
cfg.channel = label;
elec2 = ft_electrodeplacement(cfg, headshape);


%%
% demonstrate placement of 10-20 electrodes on scalp surface

fid = [
  117.9935   -2.5456   -1.5713
  11.1552  -74.6459   -8.1709
  15.4973   76.9586   -1.3693
  -78.6502    2.5375   31.1344
  ];

cfg = [];
cfg.fiducial.nas   = fid(1,:);
cfg.fiducial.ini   = fid(4,:);
cfg.fiducial.lpa   = fid(3,:);
cfg.fiducial.rpa   = fid(2,:);
cfg.method = '1020';

elec3 = ft_electrodeplacement(cfg, headshape);

%%
% demonstrate automatic placement of sEEG electrodes along shaft

% 10 electrodes placed at 5 mm distance span a total distance of 45 mm

cfg = [];
cfg.method = 'shaft';
cfg.channel = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
cfg.shaft.distance = 5;
cfg.shaft.tip   = [0 0  0]; % at the tip, i.e. deep
cfg.shaft.along = [0 0 45]; % at the end, i.e. superficial

elec0 = ft_electrodeplacement(cfg);

cfg.shaft.along = [
  8 0 22 % bend it a little bit halfway
  0 0 45 % at the end, but see below
  ];
elec1 = ft_electrodeplacement(cfg);

% note that the 10th electrodes along the curved shaft does not end up on the
% straight shaft, because the distance is compuled along the curve

figure
ft_plot_sens(elec0, 'elecshape', 'sphere', 'elecsize', 0.5, 'facecolor', 'b')
ft_plot_sens(elec1, 'elecshape', 'sphere', 'elecsize', 0.5, 'facecolor', 'r')

grid on
view(3)
cameratoolbar

%%
% demonstrate automatic placement of ECoG electrodes based on grid corners

% construct a 10x10 grid with approximately 5 mm spacing, which makes it 45x45 mm large

cfg = [];
cfg.method = 'grid';
cfg.grid.dim = [10 10];
cfg.grid.corner1 = [ 0 45 0];
cfg.grid.corner2 = [45 45 0];
cfg.grid.corner3 = [ 0  0 0];
cfg.grid.corner4 = [45  0 5]; % add some curvature by lifting corner 4 up
elec1 = ft_electrodeplacement(cfg);

elec0 = [];
elec0.elecpos(1,:) = cfg.grid.corner1;
elec0.elecpos(2,:) = cfg.grid.corner2;
elec0.elecpos(3,:) = cfg.grid.corner3;
elec0.elecpos(4,:) = cfg.grid.corner4;
elec0.label = {'1', '2', '3', '4'};

figure
ft_plot_sens(elec0, 'elecshape', 'sphere', 'elecsize', 0.7, 'facecolor', 'r')
ft_plot_sens(elec1, 'elecshape', 'sphere', 'elecsize', 0.5, 'facecolor', 'b')

grid on
view(3)
cameratoolbar

%%
% demonstrate 3D Structure Sensor headshape

headshape = ft_read_headshape(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/electrode/Model.obj'));

cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
cfg.fiducial.nas = [0.0769   -0.0064   -0.1071];
cfg.fiducial.lpa = [0.0949   -0.0119    0.0219];
cfg.fiducial.rpa = [-0.0464   -0.0070   -0.0696]; % THIS IS RATHER SLOPPY
headshape = ft_meshrealign(cfg, headshape);

cfg = [];
elec0 = ft_electrodeplacement(cfg, headshape);
