function failed_ft_prepare_headmodel

% MEM 20gb
% WALLTIME 03:00:00

% TEST test_ft_prepare_headmodel
% TEST ft_headmodel_localspheres ft_prepare_localspheres

% function to test ft_prepare_headmodel.
% This function allows for different types of structure information to be
% given, via 1) bnd, 2) seg, 3) elec, or 4) vol
% Then each of these must be dealt with appropriately for the different
% cfg.method options.
%
% (initial version by Johanna Zumer 2012)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the data which is needed

% read in the gradiometer information
hdr  = ft_read_header(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds'));
grad = hdr.grad;

% read in the mri
load standard_mri

% read in the segmented mri
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat'));

% specify the file for the headshape
shapefile  = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.shape');

% read in the headshape
shape = ft_read_headshape(shapefile);
shapevol.bnd = shape;

% hdmfile
hdmfile  = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.hdm');


%   vol = ft_prepare_headmodel(cfg)       or
%   vol = ft_prepare_headmodel(cfg, bnd)  with the output of FT_PREPARE_MESH
%   vol = ft_prepare_headmodel(cfg, seg)  with the output of FT_VOLUMESEGMENT
%   vol = ft_prepare_headmodel(cfg, elec) with the output of FT_READ_SENS
%   vol = ft_prepare_headmodel(cfg, vol)  with the output of FT_READ_VOL

csvol.o = [0, 0,0];
csvol.r = [10 50 60];
assert(ft_voltype(csvol, 'concentricspheres'))
cfg = [];
cfg.numvertices = 1000;
csbnd = ft_prepare_mesh(cfg, csvol);

ssvol.o = [0, 0,0];
ssvol.r = 60;
assert(ft_voltype(ssvol, 'singlesphere'))
cfg = [];
cfg.numvertices = 1000;
ssbnd = ft_prepare_mesh(cfg, ssvol);

cs4vol.o = [0, 0,0];
cs4vol.r = [10 30 50 60];
assert(ft_voltype(cs4vol, 'concentricspheres'))
cfg = [];
cfg.numvertices = 1000;
cs4bnd = ft_prepare_mesh(cfg, cs4vol);


load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1646/seg2'))  % brain skull scalp (3 binary masks)
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1646/seg5'))  % seg (brain, skull & scalp indexed)

elcfile  = dccnpath('/home/common/matlab/fieldtrip/data/test/original/electrodes/asa/standard_primed.elc');
elec = ft_convert_units(ft_read_sens(elcfile));


%%%%%%%%%%%%%%%%%%%%%
% do the computations

%% random/wrong cfg.method
try
  cfg = [];
  cfg.method = 'random';
  vol111 = ft_prepare_headmodel(cfg, segmentedmri);
  success = true;
catch
  success = false;
end
if success
  error('this should fail with cfg.method nonsense');
end

%% asa  - INPUT asa *.vol file

cfg = [];
cfg.hdmfile = dccnpath('/home/common/matlab/fieldtrip/data/test/original/headmodel/asa/standard.vol');
cfg.method = 'asa';
vol11 = ft_prepare_headmodel(cfg);

cfg = [];
cfg.method = 'asa';
cfg.hdmfile = hdmfile;
try % we need the asa .vol file
  vol11 = ft_prepare_headmodel(cfg);
  success = true;
catch
  success = false;
end
if success, error('11');end


%% bemcp, dipoli, openmeeg  - INPUT: segmentation or mesh (3 or also 1 compartment)

method = {'bemcp', 'dipoli','openmeeg'};
for ii = 1:numel(method)
  cfg = [];
  cfg.method = method{ii};
  
  % wrong input: nothing
  try 
    vol21 = ft_prepare_headmodel(cfg);
    success = true;
  catch
    success = false;
  end
  if success, error('21');end
  % wrong input: headshape
  try
    vol22 = ft_prepare_headmodel(cfg, shape);
    success = true;
  catch
    success = false;
  end
  if success, error('22');end
  
  % wrong input: headshape
  try
    vol23 = ft_prepare_headmodel(cfg, shapevol);
    success = true;
  catch
    success = false;
  end
  if success, error('23');end
  
  vol24 = ft_prepare_headmodel(cfg, csbnd);  % 3 compartment MESH, should work with all 3 methods
  switch method{ii}
    case 'bemcp' % this needs a 3 compartment mesh or three layered segmentation
      try % wrong input: single mesh
        vol25 = ft_prepare_headmodel(cfg, ssbnd);
        success = true;
      catch
        success = false;
      end
      if success, error('25');end
      try % wrong input: 4 compartments mesh
        vol26 = ft_prepare_headmodel(cfg, cs4bnd);
        success = true;
      catch
        success = false;
      end
      if success, error('26');end
      try % wrong input: 1 layered segmentation
        vol27 = ft_prepare_headmodel(cfg, segmentedmri);
        success = true;
      catch
        success = false;
      end
      % segmentedmari has white/gray/csf which becomes brain, but only 1 shell.
      if success, error('27');end
    case 'dipoli' % it should work also with single or 4 layers
      vol25 = ft_prepare_headmodel(cfg, ssbnd);
      vol26 = ft_prepare_headmodel(cfg, cs4bnd);
      vol27 = ft_prepare_headmodel(cfg, segmentedmri);
    case 'openmeeg' % it should work also with single compartment
      vol25 = ft_prepare_headmodel(cfg, ssbnd);
      try % but not with 4 
        vol26 = ft_prepare_headmodel(cfg, cs4bnd);
        success = true;
      catch
        success = false;
      end
      if success, error('26');end
      vol27 = ft_prepare_headmodel(cfg, segmentedmri);
  end
  
  vol28 = ft_prepare_headmodel(cfg, seg2);   % 3 layered segmentation
  vol29 = ft_prepare_headmodel(cfg, seg5);   % 3 layered segmentation (indexed style)  
  
  
end


%% concentricspheres - INPUT: 1, 3,4 compartment mesh, 1,3 layered segmentation, headshape

cfg = [];
cfg.method = 'concentricspheres';
try % wrong input: doesn't work for hdm files
  cfg.hdmfile = hdmfile;
  vol31 = ft_prepare_headmodel(cfg);
  success = true;
catch
  success = false;
end
if success, error('31');end
% it should work with everything else
% headshape
cfg = [];
cfg.method = 'concentricspheres';
cfg.headshape = shape;
vol32 = ft_prepare_headmodel(cfg);

% headshape
cfg = [];
cfg.method = 'concentricspheres';
cfg.headshape = shape.pnt;
vol33 = ft_prepare_headmodel(cfg);

% headshape
cfg = [];
cfg.method = 'concentricspheres';
cfg.headshape = shapefile;
vol34 = ft_prepare_headmodel(cfg);

% 3 comp mesh
cfg = [];
cfg.method = 'concentricspheres';
% cfg.conductivity = [.33];
vol35 = ft_prepare_headmodel(cfg, csbnd);

% 1 comp mesh
vol36 = ft_prepare_headmodel(cfg, ssbnd);
% FIXME: should we allow for this unusual but error-free way of calling
% concentricspheres?  only supplying 1 bnd and 1 conducitivity value, 
% thus output is a 1-shell vol.

% 4 comp mesh
vol37 = ft_prepare_headmodel(cfg, cs4bnd);

% 3 comp segmentation
vol38 = ft_prepare_headmodel(cfg, seg2);
vol39 = ft_prepare_headmodel(cfg, seg5);

% 1 comp segmentation
vol310 = ft_prepare_headmodel(cfg, segmentedmri);

% grad
vol311 = ft_prepare_headmodel(cfg, grad);

% elec
vol312 = ft_prepare_headmodel(cfg, elec);



%% halfspace

cfg = [];
cfg.method = 'halfspace';
% create pnt of x, y points with z on same plane (with noise)
pnt = randn(100, 3);
pnt(:, 1:2) = pnt(:,1:2)*10;
plane.pnt = pnt;
try
  vol41 = ft_prepare_headmodel(cfg, plane);
  success = true;
catch
  success = false;
end
if success, error('41');end
cfg.point = [0 0 -10];
vol42 = ft_prepare_headmodel(cfg, plane);



%% infinite

% takes no input arguments
% ft_headmodel_infinite();
cfg = [];
cfg.method = 'infinite';
vol51 = ft_prepare_headmodel(cfg);

%% localspheres
cfg         = [];
cfg.method  = 'localspheres';
try
  vol61        = ft_prepare_headmodel(cfg, shape);
  success = true;
catch
  success = false;
end
if success, error('61: grad is required for localspheres');end

cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
vol62        = ft_prepare_headmodel(cfg, shape);
vol63        = ft_prepare_headmodel(cfg, ssbnd);
try
  vol64        = ft_prepare_headmodel(cfg, csbnd);
  success = true;
catch
  success = false;
end
if success, error('64');end

vol65        = ft_prepare_headmodel(cfg, segmentedmri);
try
  vol66 = ft_prepare_headmodel(cfg, elec);
  success = true;
catch
  success = false;
end
if success, error('66');end

try % with an MRI in the input: this should fail, requires segmentation
  vol66        = ft_prepare_headmodel(cfg, mri);
  success = true;
catch
  success = false;
end
if success
  error('this should fail with unsegmented MRI input');
end
vol67 = ft_prepare_headmodel(cfg, seg2);

try % this should not work since it has 3 seg values and not specified which to use
  vol68 = ft_prepare_headmodel(cfg, seg5);
  success = true;
catch
  success = false;
end
if success, error('68');end

cfg.tissue = 3;
vol69 = ft_prepare_headmodel(cfg, seg5);


%% singleshell

cfg = [];
cfg.method = 'singleshell';
vol71 = ft_prepare_headmodel(cfg, shape);
vol72 = ft_prepare_headmodel(cfg, ssbnd);
try
  vol73 = ft_prepare_headmodel(cfg, csbnd);
  success = 1;
catch
  success = 0;
end
if success==1, error('73');end

vol74 = ft_prepare_headmodel(cfg, segmentedmri);

try
  vol75 = ft_prepare_headmodel(cfg, elec);
  success = true;
catch
  success = false;
end
if success, error('75');end

try % with an MRI in the input: this should fail, requires segmentation
  vol76        = ft_prepare_headmodel(cfg, mri);
  success = true;
catch
  success = false;
end
if success
  error('76: this should fail with unsegmented MRI input');
end

vol77 = ft_prepare_headmodel(cfg, seg2);

try % this should not work since it has 3 seg values and not specified which to use
  vol78 = ft_prepare_headmodel(cfg, seg5);
  success = true;
catch
  success = false;
end
if success, error('78');end

try
  cfg.elec = elec;
  vol79 = ft_prepare_headmodel(cfg);
  success = true;
catch
  success = false;
end
if success, error('79');end

cfg = [];
cfg.method = 'singleshell';
cfg.tissue = 3;
vol710 = ft_prepare_headmodel(cfg, seg5);


%% singlesphere

cfg = [];
cfg.method = 'singlesphere';
try
  vol81 = ft_prepare_headmodel(cfg);
  success = true;
catch
  success = false;
end
if success, error('81');end

vol82 = ft_prepare_headmodel(cfg, ssbnd);
vol83 = ft_prepare_headmodel(cfg, shape);

try
  vol84 = ft_prepare_headmodel(cfg, csbnd);
  success = 1;
catch
  success = 0;
end
if success==1, error('84');end

vol85 = ft_prepare_headmodel(cfg, segmentedmri);
vol86 = ft_prepare_headmodel(cfg, seg2);
vol87 = ft_prepare_headmodel(cfg, elec);

try
cfg.elec = elec;
vol88 = ft_prepare_headmodel(cfg);
  success = true;
catch
  success = false;
end
if success, error('88');end

try % this should not work since it has 3 seg values and not specified which to use
  vol89 = ft_prepare_headmodel(cfg, seg5);
  success = true;
catch
  success = false;
end
if success, error('89');end

cfg = [];
cfg.method = 'singlesphere';
cfg.tissue = 3;
vol810 = ft_prepare_headmodel(cfg, seg5);

cfg = [];
cfg.method = 'singlesphere';
cfg.headshape = shapefile;
vol811 = ft_prepare_headmodel(cfg);

cfg = [];
cfg.method = 'singlesphere';
cfg.headshape = shape;
vol812 = ft_prepare_headmodel(cfg);

cfg = [];
cfg.method = 'singlesphere';
cfg.headshape = shape.pnt;
vol813 = ft_prepare_headmodel(cfg);

%% simbio - INPUT: hexahedral mesh
% test_headmodel_simbio exists

clear vol*;

% you can create simbio headmodel only from hexahedral mesh
cfg = [];
cfg.method = 'hexahedral';    
hexmesh = ft_prepare_mesh(cfg, seg2);

cfg = [];
cfg.method = 'simbio';
cfg.conductivity = [0.33 0.01 0.43];  

try % this should not work 
  vol90 = ft_prepare_headmodel(cfg, mri);
  success = true;
catch
  success = false;
end
if success, error('90');end

% create headmodel
vol_hex = ft_prepare_headmodel(cfg, hexmesh);
clear vol_hex;

%% fns
% test_headmodel_fns exists

%% interpolate
% test_headmodel_interpolate

