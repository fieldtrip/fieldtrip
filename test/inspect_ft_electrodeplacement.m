function inspect_ft_electrodeplacement

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri');
mri = ft_read_mri(filename);

cfg             = [];
cfg.output      =  'scalp' ;
segmentation    = ft_volumesegment(cfg, mri);

cfg             = [];
cfg.method      = 'isosurface';
cfg.tissue      =  'scalp' ;
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

cfg = [];
cfg.method = 'volume';
cfg.channel = label;
elec1 = ft_electrodeplacement(cfg, mri);

%%

cfg = [];
cfg.method = 'headshape';
cfg.channel = label;
elec2 = ft_electrodeplacement(cfg, headshape);


%%

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
