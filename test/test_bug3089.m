function test_bug3089

% MEM 4000mb
% WALLTIME 00:20:00

% TEST test_bug3089
% TEST ft_dipolefitting ft_compute_leadfield 

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/natmeg/oddball1_mc_downsampled.fif');
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3089'));

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

cfg.continuous    = 'yes';
cfg.hpfilter      = 'no';
cfg.detrend       = 'no';
cfg.demean        = 'yes';
cfg.baselinewindow = [-inf 0];
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [50 100];
cfg.lpfilter      = 'yes';
cfg.lpfreq        = 120;
cfg.channel       = {'MEG', 'EEG'};
cfg.precision     = 'single';
cfg.coilaccuracy  = 1;
 
data_meg = ft_preprocessing(cfg);
 
%%

cfg = [];
cfg.method = 'summary';
cfg.keepchannel = 'yes';

 
cfg.channel = 'EEG';
data_meg_clean1 = ft_rejectvisual(cfg, data_meg);

 
cfg.channel = 'MEG*1';
data_meg_clean2 = ft_rejectvisual(cfg, data_meg_clean1);

 
cfg.channel = {'MEG*2', 'MEG*3'};
data_meg_clean3 = ft_rejectvisual(cfg, data_meg_clean2);

 
save data_meg_clean3 data_meg_clean3

%% reference eeg data
cfg = [];
cfg.reref = 'yes';
cfg.refchannel = 'all';
cfg.channel = 'EEG';

data_eeg = ft_preprocessing(cfg, data_meg_clean3);

cfg = [];
cfg.channel = 'MEG';

data_meg = ft_preprocessing(cfg, data_meg_clean3);

data = ft_appenddata([], data_eeg, data_meg);

%%

cfg = [];
timelock_all = ft_timelockanalysis(cfg, data);

save timelock_all timelock_all

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
timelock_cov = ft_timelockanalysis(cfg, data);

save timelock_cov timelock_cov

%% load headmodels

load headmodel_eeg
load headmodel_meg

%% convert units

headmodel_eeg = ft_convert_units(headmodel_eeg, 'm');
headmodel_meg = ft_convert_units(headmodel_meg, 'm');

timelock_all.grad  = ft_convert_units(timelock_all.grad , 'm');
timelock_all.elec  = ft_convert_units(timelock_all.elec , 'm');

timelock_cov.grad  = ft_convert_units(timelock_cov.grad , 'm');
timelock_cov.elec  = ft_convert_units(timelock_cov.elec , 'm');

%% 
cfg = [];
cfg.method = 'bemcp';
headmodel_eeg = ft_prepare_headmodel(cfg, headmodel_eeg);


%% dipole fitting 
cfg = [];
cfg.latency = 0.100;
cfg.numdipoles = 2;
cfg.symmetry = 'x';
cfg.gridsearch = 'yes';
cfg.grid.unit = 'm';
cfg.grid.resolution = 0.02;
cfg.headmodel = headmodel_meg;
cfg.channel = 'MEGMAG';
cfg.senstype = 'MEG';
dipole_mag = ft_dipolefitting(cfg, timelock_all);

cfg.headmodel = headmodel_meg;
cfg.channel = 'MEGGRAD';
dipole_grad = ft_dipolefitting(cfg, timelock_all);

cfg.headmodel = headmodel_eeg;
cfg.channel = 'EEG';
cfg.senstype = 'EEG';
dipole_eeg = ft_dipolefitting(cfg, timelock_all);

%% weighted dipole fitting
cfg = [];
cfg.channel = {'MEGMAG'};
timelock_sel = ft_selectdata(cfg, timelock_cov);

cfg.numdipoles = 2;
cfg.symmetry = 'x';
cfg.gridsearch = 'yes';
cfg.grid.unit = 'm';
cfg.grid.resolution = 0.02;
cfg.headmodel = headmodel_meg;
cfg.senstype = 'MEG';
cfg.latency = 0.100;

cfg.dipfit.noisecov = timelock_sel.cov;
dipole_mag_weighted = ft_dipolefitting(cfg, timelock_sel);


%% combined dipole fitting
cfg = [];
cfg.channel = {'MEGMAG', 'MEGGRAD'};
timelock_sel = ft_selectdata(cfg, timelock_cov);

cfg.numdipoles = 2;
cfg.symmetry = 'x';
cfg.gridsearch = 'yes';
cfg.grid.unit = 'm';
cfg.grid.resolution = 0.02;
cfg.headmodel = headmodel_meg;
cfg.senstype = 'MEG';
cfg.latency = 0.100;

cfg.dipfit.noisecov = timelock_sel.cov;
dipole_mag_grad = ft_dipolefitting(cfg, timelock_sel);

cfg.channel = {'MEGMAG', 'EEG'};
timelock_sel = ft_selectdata(cfg, timelock_cov);
dipole_mag_eeg = ft_dipolefitting(cfg, timelock_sel);


%% plot dipoles
dipoles_all = {dipole_mag, dipole_grad, dipole_eeg, dipole_mag_weighted, dipole_mag_grad};

close all
colours = {'r', 'g', 'b', 'k', 'm', 'c'};
ft_plot_axes(headmodel_meg);%, 'edgecolor', 'none', 'facealpha', 0.5);
hold on
for i = 1:length(dipoles_all)
    ft_plot_dipole(dipoles_all{i}.dip.pos, ...
        reshape(dipoles_all{i}.dip.mom, [3 2]), ...
        'unit', 'm', 'color', colours{i})
end
view([0 1 0 ])

%% plot on mri
load mri_resliced
mri_resliced = ft_convert_units(mri_resliced, 'm');
close all
figure
hold on

ft_plot_dipole(dipole_mag_weighted.dip.pos(1,:), mean(dipole_mag_weighted.dip.mom(1:3,:),2), 'color', 'r', 'unit', 'm')
ft_plot_dipole(dipole_mag_weighted.dip.pos(2,:), mean(dipole_mag_weighted.dip.mom(4:6,:),2), 'color', 'r', 'unit', 'm')

% ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g')
% ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g')

pos = mean(dipole_mag_weighted.dip.pos,1);
ft_plot_slice(mri_resliced.anatomy, 'transform',  mri_resliced.transform,  'location', pos, ...
    'orientation', [1 0 0], 'resolution', 10^0)
ft_plot_slice(mri_resliced.anatomy, transform',  mri_resliced.transform, 'location', pos, ...
    'orientation', [0 1 0], 'resolution', 10^0)
ft_plot_slice(mri_resliced.anatomy, 'transform', mri_resliced.transform, 'location', pos, ...
    'orientation', [0 0 1], 'resolution', 10^-0)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);

axis tight
axis off

view(12, -10)
