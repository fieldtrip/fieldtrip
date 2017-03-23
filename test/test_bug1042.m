function test_bug1042

% MEM 1500mb
% WALLTIME 00:45:00

% TEST test_bug1042 test_bug686
% TEST ft_convert_units ft_prepare_headmodel ft_compute_leadfield ft_prepare_sourcemodel ft_headmodel_openmeeg headsurface

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% create sensors in cm
[pnt, tri] = icosahedron162;
pnt = pnt .* 10; % convert to cm
sel = find(pnt(:,3)>0);

bnd.pnt = pnt;
bnd.unit = 'cm';

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
grad_cm = ft_convert_units(grad, 'cm');
grad_mm = ft_convert_units(grad, 'mm');

elec_m  = ft_convert_units(elec, 'm');
elec_cm = ft_convert_units(elec, 'cm');
elec_mm = ft_convert_units(elec, 'mm');

% create two conductivities
conductivity = [1 0.1];

%% For EEG the following methods are available
% singlesphere
% concentricspheres
% bem_cp
% dipoli
% openmeeg

for k = 1:numel(conductivity)
  
  cond = conductivity(k);
  
  cfg = [];
  cfg.conductivity = cond;
  
  cfg.method = 'singlesphere';
  eegvol_singlesphere = ft_prepare_headmodel(cfg, bnd);
  
  cfg.method = 'concentricspheres';
  eegvol_concentricspheres = ft_prepare_headmodel(cfg, bnd);
  % HACK otherwise it will call eeg_leadfield1 instead of eeg_leadfield4
  eegvol_concentricspheres.r    = repmat(eegvol_concentricspheres.r,    1, 4);
  eegvol_concentricspheres.cond = repmat(eegvol_concentricspheres.cond, 1, 4);
  
  bnd.pnt = pnt;
  bnd.tri = tri;
  
  geom(1) = bnd;
  geom(2) = bnd;
  geom(3) = bnd;
  geom(2).pnt = geom(2).pnt*0.9;
  geom(3).pnt = geom(3).pnt*0.8;
  
  cfg=[];
  cfg.conductivity = [cond cond cond];
  cfg.method = 'bemcp';
  eegvol_bem_cp = ft_prepare_headmodel(cfg, geom);
  
  try
    cfg.method = 'dipoli';
    eegvol_bem_dipoli = ft_prepare_headmodel(cfg, geom);
  catch
    eegvol_bem_dipoli = [];
  end
  
  try
    cfg.method = 'openmeeg';
    eegvol_bem_openmeeg = ft_prepare_headmodel(cfg, geom);
  catch
    eegvol_bem_openmeeg = [];
  end
  
  % construct them for the different geometrical units
  eegvol_singlesphere_m(k)  = ft_convert_units(eegvol_singlesphere, 'm');
  eegvol_singlesphere_cm(k) = ft_convert_units(eegvol_singlesphere, 'cm');
  eegvol_singlesphere_mm(k) = ft_convert_units(eegvol_singlesphere, 'mm');
  
  eegvol_concentricspheres_m(k)  = ft_convert_units(eegvol_concentricspheres, 'm');
  eegvol_concentricspheres_cm(k) = ft_convert_units(eegvol_concentricspheres, 'cm');
  eegvol_concentricspheres_mm(k) = ft_convert_units(eegvol_concentricspheres, 'mm');
  
  eegvol_bem_cp_m(k)  = ft_convert_units(eegvol_bem_cp, 'm');
  eegvol_bem_cp_cm(k) = ft_convert_units(eegvol_bem_cp, 'cm');
  eegvol_bem_cp_mm(k) = ft_convert_units(eegvol_bem_cp, 'mm');
  
  try
    eegvol_bem_dipoli_m(k)  = ft_convert_units(eegvol_bem_dipoli, 'm');
    eegvol_bem_dipoli_cm(k) = ft_convert_units(eegvol_bem_dipoli, 'cm');
    eegvol_bem_dipoli_mm(k) = ft_convert_units(eegvol_bem_dipoli, 'mm');
  catch
    eegvol_bem_dipoli_m(k)  = [];
    eegvol_bem_dipoli_cm(k) = [];
    eegvol_bem_dipoli_mm(k) = [];
  end
  
  try
    eegvol_bem_openmeeg_m(k)  = ft_convert_units(eegvol_bem_openmeeg, 'm');
    eegvol_bem_openmeeg_cm(k) = ft_convert_units(eegvol_bem_openmeeg, 'cm');
    eegvol_bem_openmeeg_mm(k) = ft_convert_units(eegvol_bem_openmeeg, 'mm');
  catch
    eegvol_bem_openmeeg_m(k)  = [];
    eegvol_bem_openmeeg_cm(k) = [];
    eegvol_bem_openmeeg_mm(k) = [];
  end
  
end

%% For MEG the following methods are available
% singlesphere
% localspheres
% singleshell

for k = 1:numel(conductivity)
  cond = conductivity(k);
  
  cfg = [];
  
  % FIXME this uses an undocumented option of ft_prepare_headmodel, which is due to change in the future
  cfg.conductivity = cond;
  
  cfg.method = 'singlesphere';
  megvol_singlesphere = ft_prepare_headmodel(cfg, bnd);
  
  cfg.grad = grad_cm;
  cfg.method = 'localspheres';
  megvol_localspheres = ft_prepare_headmodel(cfg, bnd);

  geom = [];
  geom.pnt = pnt;
  geom.tri = tri;
  
  cfg.method = 'singleshell';
  megvol_singleshell = ft_prepare_headmodel(cfg, geom);
  
  % construct them for the different geometrical units
  megvol_singlesphere_m(k)  = ft_convert_units(megvol_singlesphere, 'm');
  megvol_singlesphere_cm(k) = ft_convert_units(megvol_singlesphere, 'cm');
  megvol_singlesphere_mm(k) = ft_convert_units(megvol_singlesphere, 'mm');
  
  megvol_localspheres_m(k)  = ft_convert_units(megvol_localspheres, 'm');
  megvol_localspheres_cm(k) = ft_convert_units(megvol_localspheres, 'cm');
  megvol_localspheres_mm(k) = ft_convert_units(megvol_localspheres, 'mm');
  
  megvol_singleshell_m(k)  = ft_convert_units(megvol_singleshell, 'm');
  megvol_singleshell_cm(k) = ft_convert_units(megvol_singleshell, 'cm');
  megvol_singleshell_mm(k) = ft_convert_units(megvol_singleshell, 'mm');
end


%%  compute the leadfields for the volume conduction models and sensor arrays

eegvol = {
  'eegvol_singlesphere'
  'eegvol_concentricspheres'
  'eegvol_bem_cp'
  'eegvol_bem_dipoli'
  'eegvol_bem_openmeeg'
  };

megvol = {
  'megvol_singlesphere'
  'megvol_localspheres'
  'megvol_singleshell'
  };

units = {
  'm'
  'cm'
  'mm'
  };

pos = {
  [0 0 0.07]
  [0 0 7]
  [0 0 70]
  };

%%
eeg_leadfield = {};

for k = 1:numel(conductivity)
  for i=1:length(eegvol)
    for j=1:length(units)
      cfg = [];
      cfg.vol  = eval(sprintf('%s_%s(%d)', eegvol{i}, units{j}, k));
      cfg.elec = eval(sprintf('elec_%s', units{j}));
      cfg.grid.pos = pos{j};
      grid = ft_prepare_leadfield(cfg);
      eeg_leadfield{i,j,k} = grid.leadfield{1};
    end
  end
end

%%
meg_leadfield = {};
for k = 1:numel(conductivity)
  for i=1:length(megvol)
    for j=1:length(units)
      cfg      = [];
      cfg.vol  = eval(sprintf('%s_%s(%d)', megvol{i}, units{j}, k));
      cfg.grad = eval(sprintf('grad_%s', units{j}));
      cfg.grid.pos = pos{j};
      grid = ft_prepare_leadfield(cfg);
      meg_leadfield{i,j,k} = grid.leadfield{1};
    end
  end
end

%% In the table with scaling factors 
% dim1 = the "slices" correspond to different conductivities
% dim2 = the rows correspond to the different volume conduction models
% dim3 = the columns correspond to m, cm, mm

eeg_table = cellfun(@norm, eeg_leadfield);
meg_table = cellfun(@norm, meg_leadfield);

format long
(eeg_table ./ eeg_table(1,1));
(meg_table ./ meg_table(1,1));
format short

