function test_bug3089

% MEM 4000mb
% WALLTIME 00:20:00

% TEST ft_dipolefitting ft_compute_leadfield ft_apply_transform

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/natmeg/oddball1_mc_downsampled.fif');
datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3089');
addpath(datadir); % for the trialfuns

%%

cfg = [];
cfg.dataset = dataset;
cfg.trialdef.prestim        = 0.2;
cfg.trialdef.poststim       = 0.4;
cfg.trialdef.rsp_triggers   = [256 4096];
cfg.trialdef.stim_triggers  = [1 2];
cfg.trialfun                = 'trialfun_oddball_stimlocked';

cfg = ft_definetrial(cfg);

%%

cfg.continuous     = 'yes';
cfg.hpfilter       = 'no';
cfg.detrend        = 'no';
cfg.demean         = 'yes';
cfg.baselinewindow = [-inf 0];
cfg.dftfilter      = 'yes';
cfg.dftfreq        = [50 100];
cfg.lpfilter       = 'yes';
cfg.lpfreq         = 120;
cfg.channel        = {'MEG', 'EEG'};
cfg.precision      = 'single';
cfg.coilaccuracy   = 1; % use the COILDEF file

data_raw = ft_preprocessing(cfg);


%% reject noisy trials

if exist(fullfile(datadir,'data_clean.mat'), 'file')
  load(fullfile(datadir, 'data_clean.mat'));

else
  cfg = [];
  cfg.method = 'summary';
  cfg.keepchannel = 'yes';
  
  cfg.channel = 'EEG';
  data_clean = ft_rejectvisual(cfg, data_raw);
  
  cfg.channel = 'MEG*1'; % MEGMAG
  data_clean = ft_rejectvisual(cfg, data_clean);
  
  cfg.channel = {'MEG*2', 'MEG*3'}; % MEGGRAD
  data_clean = ft_rejectvisual(cfg, data_clean);
  
  save(fullfile(datadir,'data_clean'), 'data_clean');
end

%% reference eeg data

if false
  cfg = [];
  cfg.reref = 'yes';
  cfg.refchannel = 'all';
  cfg.channel = 'EEG';
  data_eeg = ft_preprocessing(cfg, data_clean);
  
  cfg = [];
  cfg.channel = 'MEG';
  data_meg = ft_preprocessing(cfg, data_clean);
  
  cfg = []
  data_all = ft_appenddata(cfg, data_eeg, data_meg);
  
else
  montage = [];
  montage.labelold = ft_channelselection('EEG', data_clean.label);
  montage.labelnew = ft_channelselection('EEG', data_clean.label);
  montage.tra = eye(length(montage.labelnew), length(montage.labelold));
  for i=1:length(montage.labelnew)
    montage.tra(i,:) = montage.tra(i,:) - ones(1,length(montage.labelold))/length(montage.labelold);
  end
  data_all      = ft_apply_montage(data_clean, montage, 'keepunused', true, 'balancename', 'avgref');
  
  % apply the montage to the electrode definition
  elec = data_clean.elec;
  elec.tra = eye(length(elec.label));
  elec.balance.current = 'none';
  elec = ft_apply_montage(elec, montage, 'keepunused', true, 'balancename', 'avgref');
  
  % keep it with the data
  data_all.elec = elec;
end


%% average

cfg = [];
timelock_all = ft_timelockanalysis(cfg, data_all);

save(fullfile(datadir,'timelock_all'),'timelock_all');

%%

cfg = [];
cfg.layout = 'neuromag306mag.lay';
figure; ft_multiplotER(cfg, timelock_all);

cfg.layout = 'neuromag306planar.lay';
figure; ft_multiplotER(cfg, timelock_all);

cfg = [];
cfg.elec = timelock_all.elec;
cfg.layout = ft_prepare_layout(cfg);
figure; ft_multiplotER(cfg, timelock_all);

%%

cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-0.2 0];
timelock_cov = ft_timelockanalysis(cfg, data_all);

save(fullfile(datadir,'timelock_cov'),'timelock_cov');

%% load headmodels

load(fullfile(datadir,'headmodel_eeg'));
load(fullfile(datadir,'headmodel_meg'));

%% convert units

headmodel_eeg = ft_convert_units(headmodel_eeg, 'm');
headmodel_meg = ft_convert_units(headmodel_meg, 'm');

timelock_all.grad  = ft_convert_units(timelock_all.grad , 'm');
timelock_all.elec  = ft_convert_units(timelock_all.elec , 'm');

timelock_cov.grad  = ft_convert_units(timelock_cov.grad , 'm');
timelock_cov.elec  = ft_convert_units(timelock_cov.elec , 'm');

%% recompute the EEG headmodel (following unit change)

cfg = [];
cfg.method = 'bemcp';
headmodel_eeg = ft_prepare_headmodel(cfg, headmodel_eeg);

%% conventional dipole fitting

cfg = [];
cfg.latency         = 0.100;
cfg.numdipoles      = 2;
cfg.symmetry        = 'x';
cfg.gridsearch      = 'yes';
cfg.grid.unit       = 'm';
cfg.grid.resolution = 0.02;

cfg.senstype        = 'MEG';
cfg.headmodel       = headmodel_meg;
cfg.channel         = 'MEGMAG';
dipole_mag = ft_dipolefitting(cfg, timelock_all);

cfg.headmodel       = headmodel_meg;
cfg.channel         = 'MEGGRAD';
dipole_grad = ft_dipolefitting(cfg, timelock_all);

cfg.senstype        = 'EEG';
cfg.headmodel       = headmodel_eeg;
cfg.channel         = 'EEG';
dipole_eeg = ft_dipolefitting(cfg, timelock_all);

%% weighted dipole fitting

cfg = [];

cfg.latency         = 0.100;
cfg.numdipoles      = 2;
cfg.symmetry        = 'x';
cfg.gridsearch      = 'yes';
cfg.grid.unit       = 'm';
cfg.grid.resolution = 0.02;

cfg.channel = {'MEGMAG'};
timelock_sel = ft_selectdata(cfg, timelock_cov); % this selects channels from the covariance
cfg.senstype        = 'MEG';
cfg.headmodel       = headmodel_meg;
cfg.dipfit.noisecov = timelock_sel.cov;
dipole_mag_weighted = ft_dipolefitting(cfg, timelock_sel);

cfg.channel = {'MEGGRAD'};
timelock_sel = ft_selectdata(cfg, timelock_cov); % this selects channels from the covariance
cfg.senstype        = 'MEG';
cfg.headmodel       = headmodel_meg;
cfg.dipfit.noisecov = timelock_sel.cov;
dipole_grad_weighted = ft_dipolefitting(cfg, timelock_sel);

cfg.channel = {'EEG'};
timelock_sel = ft_selectdata(cfg, timelock_cov); % this selects channels from the covariance
cfg.senstype        = 'EEG';
cfg.headmodel       = headmodel_eeg;
cfg.dipfit.noisecov = timelock_sel.cov;
dipole_eeg_weighted = ft_dipolefitting(cfg, timelock_sel);


%% weighted combined dipole fitting

cfg = [];

cfg.numdipoles      = 2;
cfg.latency         = 0.100;
cfg.symmetry        = 'x';
cfg.gridsearch      = 'yes';
cfg.grid.unit       = 'm';
cfg.grid.resolution = 0.02;

cov_mag  = zeros(306);
for i=1:3:306
  cov_mag(i,i) = 1;
end

cov_grad = zeros(306);
for i=2:3:306
  cov_grad(i,i) = 1;
end
for i=3:3:306
  cov_grad(i,i) = 1;
end


cov_grad = zeros(306);
for i=2:3:306
  cov_grad(i,i) = 1;
end
for i=3:3:306
  cov_grad(i,i) = 1;
end

cov_58 = zeros(306);
for i=1:3:306
  cov_58(i,i) = 1;
end
for i=2:3:306
  cov_58(i,i) = 58^2;
end
for i=3:3:306
  cov_58(i,i) = 58^2;
end

cov_20 = zeros(306);
for i=1:3:306
  cov_20(i,i) = 1;
end
for i=2:3:306
  cov_20(i,i) = 20^2;
end
for i=3:3:306
  cov_20(i,i) = 20^2;
end


cfg.channel = {'MEGMAG', 'MEGGRAD'};
timelock_sel = ft_selectdata(cfg, timelock_cov); % this selects channels from the covariance
cfg.senstype        = 'MEG';
cfg.headmodel       = headmodel_meg;
cfg.dipfit.noisecov = timelock_sel.cov;
% cfg.dipfit.noisecov = cov_20
dipole_mag_grad = ft_dipolefitting(cfg, timelock_sel);

cfg.channel = {'MEGMAG', 'EEG'};
timelock_sel = ft_selectdata(cfg, timelock_cov);
cfg.senstype        = {'MEG', 'EEG'};
cfg.headmodel       = {headmodel_meg, headmodel_eeg};
cfg.dipfit.noisecov = timelock_sel.cov;
dipole_mag_eeg = ft_dipolefitting(cfg, timelock_sel);

cfg.channel = {'MEGGRAD', 'EEG'};
timelock_sel = ft_selectdata(cfg, timelock_cov);
cfg.senstype        = {'MEG', 'EEG'};
cfg.headmodel       = {headmodel_meg, headmodel_eeg};
cfg.dipfit.noisecov = timelock_sel.cov;
dipole_grad_eeg = ft_dipolefitting(cfg, timelock_sel);

cfg.channel = {'MEGMAG', 'MEGGRAD', 'EEG'};
timelock_sel = ft_selectdata(cfg, timelock_cov);
cfg.senstype        = {'MEG', 'EEG'};
cfg.headmodel       = {headmodel_meg, headmodel_eeg};
cfg.dipfit.noisecov = timelock_sel.cov;
dipole_mag_grad_eeg = ft_dipolefitting(cfg, timelock_sel);


%% plot dipoles

dipoles_all = {dipole_mag, dipole_grad, dipole_eeg, dipole_mag_grad};

close all
colours = {'r', 'g', 'b', 'k', 'm', 'c'};
ft_plot_axes(headmodel_meg);%, 'edgecolor', 'none', 'facealpha', 0.5);
hold on
for i = 1:length(dipoles_all)
  ft_plot_dipole(dipoles_all{i}.dip.pos, reshape(dipoles_all{i}.dip.mom, [3 2]), 'unit', 'm', 'color', colours{i})
end
view([0 1 0 ])

%% plot dipoles on mri

load mri_resliced
mri_resliced = ft_convert_units(mri_resliced, 'm');

close all
figure
hold on

for i = 1:length(dipoles_all)
  ft_plot_dipole(dipoles_all{i}.dip.pos, reshape(dipoles_all{i}.dip.mom, [3 2]), 'unit', 'm', 'color', colours{i})
end

pos = mean(dipole_mag_weighted.dip.pos,1);
ft_plot_slice(mri_resliced.anatomy, 'transform', mri_resliced.transform, 'location', pos, 'orientation', [1 0 0]);
ft_plot_slice(mri_resliced.anatomy, 'transform', mri_resliced.transform, 'location', pos, 'orientation', [0 1 0]);
ft_plot_slice(mri_resliced.anatomy, 'transform', mri_resliced.transform, 'location', pos, 'orientation', [0 0 1]);

ft_plot_crosshair(pos, 'color', [1 1 1]/2);

axis auto
axis tight
axis off

view(12, -10)
