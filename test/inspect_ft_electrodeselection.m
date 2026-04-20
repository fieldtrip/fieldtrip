function inspect_ft_electrodeselection

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_electrodeplacement ft_electrodeselection

cd(dccnpath('/project/3031000.02/external/download/tutorial/opm_helmet_design'))

%%

headshape = ft_read_headshape('spherical-head.stl');
helmet = ft_read_headshape('spherical-helmet.stl');

% the CTF coordinate system has X towards the nose and Y towards the left
headshape.coordsys = 'ctf';
helmet.coordsys = 'ctf';

%%

nas = [+100 0 0];
ini = [-100 0 0];
lpa = [0 +100 0];
rpa = [0 -100 0];

headshape.fid.pos = [
  nas
  ini
  lpa
  rpa
  ];

headshape.fid.label = {
  'nas'
  'ini'
  'lpa'
  'rpa'
  };

%%

% this places a lot of electrodes, see https://doi.org/10.1016/s1388-2457(00)00527-7
cfg = [];
cfg.fiducial.nas   = nas;
cfg.fiducial.ini   = ini;
cfg.fiducial.lpa   = lpa;
cfg.fiducial.rpa   = rpa;
cfg.method = '1020';
cfg.feedback = 'no';
elec = ft_electrodeplacement(cfg, headshape);

%%

cfg = [];
cfg.elec = elec;
cfg.outwardshift = 10/2; % the helmet is 10 mm thick, the bottom of the sensor holder will be halfway in

cfg.template = 'fieldline_holder.stl';
[tmpcfg, holder] = ft_sensorplacement(cfg, headshape);

%%

% this makes a selection only on the basis of the labels

cfg = [];
cfg.channel = 'eeg1020';
elec1020 = ft_electrodeselection(cfg, elec);

%%

% this makes a selection with the graphical user interface

cfg = [];
cfg.headshape = headshape;
cfg.mesh = holder;
cfg.headshapestyle.material = 'dull'; % see FT_PLOT_HEADSHAPE
cfg.meshstyle = [];                   % see FT_PLOT_MESH
cfg.elecstyle.elecshape = 'disc';     % see FT_PLOT_SENS
cfg.channel = elec1020.label;         % the initial selection
selected = ft_electrodeselection(cfg, elec);
