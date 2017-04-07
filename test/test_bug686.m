function test_bug686

% MEM 2gb
% WALLTIME 00:30:00

% TEST ft_convert_units ft_prepare_headmodel ft_prepare_leadfield ft_prepare_sourcemodel ft_headmodel_openmeeg headsurface

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

[pnt, tri] = icosahedron162;
pnt = pnt .* 10;           % convert to cm
sel = find(pnt(:,3)>0);    % take the upper hemisphere

elec.pnt = pnt(sel,:);
for i=1:length(sel)
  elec.label{i} = sprintf('electrode%d', i);
end
elec.unit = 'cm';

grad.pnt = pnt(sel,:) .* 1.2;
grad.ori = pnt(sel,:);
grad.tra = eye(length(sel));
for i=1:length(sel)
  grad.ori(i,:) = grad.ori(i,:) ./ norm(grad.ori(i,:));
  grad.label{i} = sprintf('magnetometer%d', i);
end
grad.unit = 'cm';

% construct them for the different geometrical units
grad_m  = ft_convert_units(grad, 'm');
grad_dm = ft_convert_units(grad, 'dm');
grad_cm = ft_convert_units(grad, 'cm');
grad_mm = ft_convert_units(grad, 'mm');

elec_m  = ft_convert_units(elec, 'm');
elec_dm = ft_convert_units(elec, 'dm');
elec_cm = ft_convert_units(elec, 'cm');
elec_mm = ft_convert_units(elec, 'mm');


%% For EEG the following methods are available

% for some of them I cannot test with a sphere
%  cfg.method = 'bem_asa'
%  cfg.method = 'halfspace
%  cfg.method = 'infinite

geom1 = [];
geom1.pnt = pnt;
geom1.tri = tri;

geom3 = [];
geom3(1).tri = tri;
geom3(2).tri = tri;
geom3(3).tri = tri;
geom3(1).pnt = pnt;
geom3(2).pnt = pnt*0.9;
geom3(3).pnt = pnt*0.8;

cfg = [];
cfg.conductivity = 1;
cfg.method = 'singlesphere';
eegvol_singlesphere = ft_prepare_headmodel(cfg, geom1);
eegvol_singlesphere.o = [0 0 0]; % avoid rounding off errors

cfg.conductivity = [1 1 1];
cfg.method = 'concentricspheres';
eegvol_concentricspheres = ft_prepare_headmodel(cfg, geom3);
eegvol_concentricspheres.o = [0 0 0]; % avoid rounding off errors

cfg.method = 'bemcp'; % this is only implemented for 3 compartments
cfg.conductivity = [1 1 1];
eegvol_bemcp3 = ft_prepare_headmodel(cfg, geom3);

% some of the fwd solutions require the external toolbox, which is not available on all platforms
try
  cfg.method = 'dipoli';
  cfg.conductivity = 1;
  eegvol_dipoli1 = ft_prepare_headmodel(cfg, geom1);
  cfg.conductivity = [1 1 1];
  eegvol_dipoli3 = ft_prepare_headmodel(cfg, geom3);
catch
  eegvol_dipoli1 = [];
  eegvol_dipoli3 = [];
end

try 
  cfg.method = 'openmeeg';
  cfg.conductivity = 1;
  eegvol_openmeeg1 = ft_prepare_headmodel(cfg, geom1);
  cfg.conductivity = [1 1 1];
  eegvol_openmeeg3 = ft_prepare_headmodel(cfg, geom3);
catch
  eegvol_openmeeg1 = [];
  eegvol_openmeeg3 = [];
end

% construct them for the different geometrical units
eegvol_singlesphere_m  = ft_convert_units(eegvol_singlesphere, 'm');
eegvol_singlesphere_dm = ft_convert_units(eegvol_singlesphere, 'dm');
eegvol_singlesphere_cm = ft_convert_units(eegvol_singlesphere, 'cm');
eegvol_singlesphere_mm = ft_convert_units(eegvol_singlesphere, 'mm');

eegvol_concentricspheres_m  = ft_convert_units(eegvol_concentricspheres, 'm');
eegvol_concentricspheres_dm = ft_convert_units(eegvol_concentricspheres, 'dm');
eegvol_concentricspheres_cm = ft_convert_units(eegvol_concentricspheres, 'cm');
eegvol_concentricspheres_mm = ft_convert_units(eegvol_concentricspheres, 'mm');

eegvol_bemcp3_m  = ft_convert_units(eegvol_bemcp3, 'm');
eegvol_bemcp3_dm = ft_convert_units(eegvol_bemcp3, 'dm');
eegvol_bemcp3_cm = ft_convert_units(eegvol_bemcp3, 'cm');
eegvol_bemcp3_mm = ft_convert_units(eegvol_bemcp3, 'mm');

try
  eegvol_dipoli1_m  = ft_convert_units(eegvol_dipoli1, 'm');
  eegvol_dipoli1_dm = ft_convert_units(eegvol_dipoli1, 'dm');
  eegvol_dipoli1_cm = ft_convert_units(eegvol_dipoli1, 'cm');
  eegvol_dipoli1_mm = ft_convert_units(eegvol_dipoli1, 'mm');
  eegvol_dipoli3_m  = ft_convert_units(eegvol_dipoli3, 'm');
  eegvol_dipoli3_dm = ft_convert_units(eegvol_dipoli3, 'dm');
  eegvol_dipoli3_cm = ft_convert_units(eegvol_dipoli3, 'cm');
  eegvol_dipoli3_mm = ft_convert_units(eegvol_dipoli3, 'mm');
catch
  fprintf('Please install Dipoli\n')
  % leaving them empty will be interpreted as an infinite volume conductor
  eegvol_dipoli1_m  = [];
  eegvol_dipoli1_dm = [];
  eegvol_dipoli1_cm = [];
  eegvol_dipoli1_mm = [];  
  eegvol_dipoli3_m  = [];
  eegvol_dipoli3_dm = [];
  eegvol_dipoli3_cm = [];
  eegvol_dipoli3_mm = [];  
end

try 
  eegvol_openmeeg1_m  = ft_convert_units(eegvol_openmeeg1, 'm');
  eegvol_openmeeg1_dm = ft_convert_units(eegvol_openmeeg1, 'dm');
  eegvol_openmeeg1_cm = ft_convert_units(eegvol_openmeeg1, 'cm');
  eegvol_openmeeg1_mm = ft_convert_units(eegvol_openmeeg1, 'mm');  
  eegvol_openmeeg3_m  = ft_convert_units(eegvol_openmeeg3, 'm');
  eegvol_openmeeg3_dm = ft_convert_units(eegvol_openmeeg3, 'dm');
  eegvol_openmeeg3_cm = ft_convert_units(eegvol_openmeeg3, 'cm');
  eegvol_openmeeg3_mm = ft_convert_units(eegvol_openmeeg3, 'mm');  
catch
  fprintf('Please install OpenMEEG\n')
  % leaving them empty will be interpreted as an infinite volume conductor
  eegvol_openmeeg1_m  = [];
  eegvol_openmeeg1_dm = [];
  eegvol_openmeeg1_cm = [];
  eegvol_openmeeg1_mm = [];
  eegvol_openmeeg3_m  = [];
  eegvol_openmeeg3_dm = [];
  eegvol_openmeeg3_cm = [];
  eegvol_openmeeg3_mm = [];
end

%% For MEG the following methods are available

% for some of them I cannot test with a sphere
% cfg.method = 'infinite'

geom = [];
geom.pnt = pnt;

cfg = [];
cfg.conductivity = 1;

cfg.method = 'singlesphere';
megvol_singlesphere = ft_prepare_headmodel(cfg, geom);

cfg.grad = grad_cm;
cfg.method = 'localspheres';
megvol_localspheres = ft_prepare_headmodel(cfg, geom);

geom.tri = tri;
cfg.method = 'singleshell';
megvol_singleshell = ft_prepare_headmodel(cfg, geom);

% construct them for the different geometrical units
megvol_singlesphere_m  = ft_convert_units(megvol_singlesphere, 'm');
megvol_singlesphere_dm = ft_convert_units(megvol_singlesphere, 'dm');
megvol_singlesphere_cm = ft_convert_units(megvol_singlesphere, 'cm');
megvol_singlesphere_mm = ft_convert_units(megvol_singlesphere, 'mm');

megvol_localspheres_m  = ft_convert_units(megvol_localspheres, 'm');
megvol_localspheres_dm = ft_convert_units(megvol_localspheres, 'dm');
megvol_localspheres_cm = ft_convert_units(megvol_localspheres, 'cm');
megvol_localspheres_mm = ft_convert_units(megvol_localspheres, 'mm');

megvol_singleshell_m  = ft_convert_units(megvol_singleshell, 'm');
megvol_singleshell_dm = ft_convert_units(megvol_singleshell, 'dm');
megvol_singleshell_cm = ft_convert_units(megvol_singleshell, 'cm');
megvol_singleshell_mm = ft_convert_units(megvol_singleshell, 'mm');


%%  compute the leadfields for the volume conduction models and sensor arrays

eegvol = {
  'eegvol_singlesphere'
  'eegvol_concentricspheres'
  'eegvol_bemcp3'
  'eegvol_dipoli1'
  'eegvol_dipoli3'
  'eegvol_openmeeg1'
  'eegvol_openmeeg3'
  };

megvol = {
  'megvol_singlesphere'
  'megvol_localspheres'
  'megvol_singleshell'
  };

units = {
  'm'
  'dm'
  'cm'
  'mm'
  };

pos = {
  [0 0  0.07]
  [0 0  0.70]
  [0 0  7.00]
  [0 0 70.00]
  };

%%
eeg_leadfield = {};
for i=1:length(eegvol)
  for j=1:length(units)
    cfg             = [];
    cfg.vol         = eval(sprintf('%s_%s', eegvol{i}, units{j}));
    cfg.elec        = eval(sprintf('elec_%s', units{j}));
    cfg.grid.pos    = pos{j};
    cfg.sourceunits = units{j};
    grid = ft_prepare_leadfield(cfg);
    eeg_leadfield{i,j} = grid.leadfield{1};
  end
end

%%
meg_leadfield = {};
for i=1:length(megvol)
  for j=1:length(units)
    cfg             = [];
    cfg.vol         = eval(sprintf('%s_%s', megvol{i}, units{j}));
    cfg.grad        = eval(sprintf('grad_%s', units{j}));
    cfg.grid.pos    = pos{j};
    cfg.sourceunits = units{j};
    grid = ft_prepare_leadfield(cfg);
    meg_leadfield{i,j} = grid.leadfield{1};
  end
end

%% In the table with scaling factors the columns correspond to m, cm, mm, 
% the rows correspond to the different volume conduction models

eeg_table = cellfun(@norm, eeg_leadfield);
disp('eeg_table')
disp(round(log10(eeg_table ./ eeg_table(1,1))))

disp('meg_table')
meg_table = cellfun(@norm, meg_leadfield);
disp(round(log10(meg_table ./ meg_table(1,1))))


