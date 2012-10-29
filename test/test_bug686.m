function test_bug686

% TEST test_bug686
% TEST ft_convert_units ft_prepare_headmodel ft_prepare_leadfield

global ft_default
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
grad_cm = ft_convert_units(grad, 'cm');
grad_mm = ft_convert_units(grad, 'mm');

elec_m  = ft_convert_units(elec, 'm');
elec_cm = ft_convert_units(elec, 'cm');
elec_mm = ft_convert_units(elec, 'mm');


%% For EEG the following methods are available

% for some of them I cannot test with a sphere
%  cfg.method = 'bem_asa'
%  cfg.method = 'halfspace
%  cfg.method = 'infinite

cfg = [];

% FIXME this uses an undocumented option of ft_prepare_headmodel, which is due to change in the future
cfg.geom.pnt = pnt;
cfg.conductivity = 1;

cfg.method = 'singlesphere';
eegvol_singlesphere = ft_prepare_headmodel(cfg);

cfg.method = 'concentricspheres';
eegvol_concentricspheres = ft_prepare_headmodel(cfg);
% HACK otherwise it will call eeg_leadfield1 instead of eeg_leadfield4
eegvol_concentricspheres.r = repmat(eegvol_concentricspheres.r, 1, 4);
eegvol_concentricspheres.c = repmat(eegvol_concentricspheres.c, 1, 4);


bnd.pnt = pnt;
bnd.tri = tri;

cfg=[];
cfg.geom(1) = bnd;
cfg.geom(2) = bnd;
cfg.geom(3) = bnd;
cfg.geom(2).pnt = cfg.geom(2).pnt*0.9;
cfg.geom(3).pnt = cfg.geom(3).pnt*0.8;
cfg.conductivity = [1 1 1];
cfg.method = 'bem_cp';
eegvol_bem_cp = ft_prepare_headmodel(cfg);

% some of the fwd solutions require the external toolbox, not always
% downloaded
try
  cfg.method = 'bem_dipoli';
  eegvol_bem_dipoli = ft_prepare_headmodel(cfg);
catch
  eegvol_bem_dipoli = [];
end

try 
  cfg.method = 'bem_openmeeg';
  eegvol_bem_openmeeg = ft_prepare_headmodel(cfg);
catch
  eegvol_bem_openmeeg = [];
end

% construct them for the different geometrical units
eegvol_singlesphere_m  = ft_convert_units(eegvol_singlesphere, 'm');
eegvol_singlesphere_cm = ft_convert_units(eegvol_singlesphere, 'cm');
eegvol_singlesphere_mm = ft_convert_units(eegvol_singlesphere, 'mm');

eegvol_concentricspheres_m  = ft_convert_units(eegvol_concentricspheres, 'm');
eegvol_concentricspheres_cm = ft_convert_units(eegvol_concentricspheres, 'cm');
eegvol_concentricspheres_mm = ft_convert_units(eegvol_concentricspheres, 'mm');

eegvol_bem_cp_m  = ft_convert_units(eegvol_bem_cp, 'm');
eegvol_bem_cp_cm = ft_convert_units(eegvol_bem_cp, 'cm');
eegvol_bem_cp_mm = ft_convert_units(eegvol_bem_cp, 'mm');

try
  eegvol_bem_dipoli_m  = ft_convert_units(eegvol_bem_dipoli, 'm');
  eegvol_bem_dipoli_cm = ft_convert_units(eegvol_bem_dipoli, 'cm');
  eegvol_bem_dipoli_mm = ft_convert_units(eegvol_bem_dipoli, 'mm');
catch
  fprintf('Please install Dipoli')
  % leaving them empty will be interpreted as an infinite volume conductor
  eegvol_bem_dipoli_m  = [];
  eegvol_bem_dipoli_cm = [];
  eegvol_bem_dipoli_mm = [];  
end

try 
  eegvol_bem_openmeeg_m  = ft_convert_units(eegvol_bem_openmeeg, 'm');
  eegvol_bem_openmeeg_cm = ft_convert_units(eegvol_bem_openmeeg, 'cm');
  eegvol_bem_openmeeg_mm = ft_convert_units(eegvol_bem_openmeeg, 'mm');  
catch
  fprintf('Please install OpenMEEG')
  % leaving them empty will be interpreted as an infinite volume conductor
  eegvol_bem_openmeeg_m  = [];
  eegvol_bem_openmeeg_cm = [];
  eegvol_bem_openmeeg_mm = [];
end

%% For MEG the following methods are available

% for some of them I cannot test with a sphere
% cfg.method = 'infinite'

cfg = [];

% FIXME this uses an undocumented option of ft_prepare_headmodel, which is due to change in the future
cfg.geom.pnt = pnt;
cfg.conductivity = 1;

cfg.method = 'singlesphere';
megvol_singlesphere = ft_prepare_headmodel(cfg);

cfg.grad = grad_cm;
cfg.method = 'localspheres';
megvol_localspheres = ft_prepare_headmodel(cfg);

cfg.geom.tri = tri;
cfg.method = 'singleshell';
megvol_singleshell = ft_prepare_headmodel(cfg);

% construct them for the different geometrical units
megvol_singlesphere_m  = ft_convert_units(megvol_singlesphere, 'm');
megvol_singlesphere_cm = ft_convert_units(megvol_singlesphere, 'cm');
megvol_singlesphere_mm = ft_convert_units(megvol_singlesphere, 'mm');

megvol_localspheres_m  = ft_convert_units(megvol_localspheres, 'm');
megvol_localspheres_cm = ft_convert_units(megvol_localspheres, 'cm');
megvol_localspheres_mm = ft_convert_units(megvol_localspheres, 'mm');

megvol_singleshell_m  = ft_convert_units(megvol_singleshell, 'm');
megvol_singleshell_cm = ft_convert_units(megvol_singleshell, 'cm');
megvol_singleshell_mm = ft_convert_units(megvol_singleshell, 'mm');


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
for i=1:length(eegvol)
  for j=1:length(units)
    cfg = [];
    cfg.vol  = eval(sprintf('%s_%s', eegvol{i}, units{j}));
    cfg.elec = eval(sprintf('elec_%s', units{j}));
    cfg.grid.pos = pos{j};
    grid = ft_prepare_leadfield(cfg);
    eeg_leadfield{i,j} = grid.leadfield{1};
  end
end

%%
meg_leadfield = {};
for i=1:length(megvol)
  for j=1:length(units)
    cfg      = [];
    cfg.vol  = eval(sprintf('%s_%s', megvol{i}, units{j}));
    cfg.grad = eval(sprintf('grad_%s', units{j}));
    cfg.grid.pos = pos{j};
    grid = ft_prepare_leadfield(cfg);
    meg_leadfield{i,j} = grid.leadfield{1};
  end
end

%% In the table with scaling factors the columns correspond to m, cm, mm, 
% the rows correspond to the different volume conduction models

eeg_table = cellfun(@norm, eeg_leadfield);
meg_table = cellfun(@norm, meg_leadfield);

log10(eeg_table ./ eeg_table(1,1))
log10(meg_table ./ meg_table(1,1))

return
% FIXME the rest is not used
% 
% Loop over all the methods and express the median ratio of first/second and
% second/third units leadfields, for all components

%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG
%%%%%%%%%%%%%%%%%%%%%%%%%

for eeg_method =1:size(eeg_leadfield,1)
  eeg_method
  ori = 1;
  sprintf('Component 1')
  tt=[eeg_leadfield{eeg_method,1}(:,ori) eeg_leadfield{eeg_method,2}(:,ori) eeg_leadfield{eeg_method,3}(:,ori)];
  disp(['first/second unit median ratio for EEG:'])
  median(tt(:,1)./tt(:,2))
  disp(['second/third unit median ratio for EEG:'])
  median(tt(:,2)./tt(:,3))
  
  ori = 2;
  sprintf('Component 2')
  tt=[eeg_leadfield{eeg_method,1}(:,ori) eeg_leadfield{eeg_method,2}(:,ori) eeg_leadfield{eeg_method,3}(:,ori)];
  disp(['first/second unit median ratio for EEG:'])
  median(tt(:,1)./tt(:,2))
  disp(['second/third unit median ratio for EEG:'])
  median(tt(:,2)./tt(:,3))
  
  ori = 3;
  sprintf('Component 3')
  tt=[eeg_leadfield{eeg_method,1}(:,ori) eeg_leadfield{eeg_method,2}(:,ori) eeg_leadfield{eeg_method,3}(:,ori)];
  disp(['first/second unit median ratio for EEG:'])
  median(tt(:,1)./tt(:,2))
  disp(['second/third unit median ratio for EEG:'])
  median(tt(:,2)./tt(:,3))
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% MEG
%%%%%%%%%%%%%%%%%%%%%%%%%

for meg_method =1:size(meg_leadfield,1)
  meg_method
  % the third component for MEG is not reliable, so we test only ori=1 and ori=2
  ori = 1; 
  sprintf('Component 1')
  tt=[meg_leadfield{meg_method,1}(:,ori) meg_leadfield{meg_method,2}(:,ori) meg_leadfield{meg_method,3}(:,ori)];
  disp(['first/second unit median ratio for EEG:'])
  median(tt(:,1)./tt(:,2))
  disp(['second/third unit median ratio for EEG:'])
  median(tt(:,2)./tt(:,3))

  ori = 2; 
  sprintf('Component 2')
  tt=[meg_leadfield{meg_method,1}(:,ori) meg_leadfield{meg_method,2}(:,ori) meg_leadfield{meg_method,3}(:,ori)];
  disp(['first/second unit median ratio for EEG:'])
  median(tt(:,1)./tt(:,2))
  disp(['second/third unit median ratio for EEG:'])
  median(tt(:,2)./tt(:,3))
  
end
