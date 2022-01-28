function test_example_compute_forward_simulated_data

% MEM 4gb
% WALLTIME 00:10:00

%
%% Compute forward simulated data using ft_dipolesimulation
%
%%
% load the template MNI headmodel, electrodes and MRI

[ftver, ftpath] = ft_version;
headmodel = ft_read_headmodel(fullfile(ftpath, 'template/headmodel/standard_bem.mat'));
elec = ft_read_sens(fullfile(ftpath, 'template/electrode/standard_1020.elc'));
mri = ft_read_mri(fullfile(ftpath, 'template/anatomy/single_subj_T1_1mm.nii'));

%%
% explore the location to place the dipole

cfg = [];
ft_sourceplot(cfg, mri);

%
%%
% compute a forward model for a single dipole, 10 trials of 1 second each, with a 2Hz signal

cfg = [];
cfg.dip.unit = 'mm';
cfg.dip.pos = [-40 -20 50]; % left motor cortex
% cfg.dip.pos = [-50 30 -10]; % orbitofrontaal
cfg.dip.mom = cfg.dip.pos/norm(cfg.dip.pos); % radial
% cfg.dip.mom = [0 1 0]; % tangential
cfg.dip.frequency = 2;
cfg.elec = elec;
cfg.headmodel = headmodel;
data = ft_dipolesimulation(cfg);

%%
% use low-level functions to make a detailled figure

figure
% ft_plot_headmodel(headmodel);
% ft_plot_ortho(mri.anatomy, 'location', [0 0 0], 'transform', mri.tra);
ft_plot_sens(elec, 'label', 'label');
ft_plot_dipole(data.cfg.dip.pos, data.cfg.dip.mom, 'unit', 'mm')

%%
% again look at ft_sourceplot, make a cross-section at the dipole position

cfg = [];
cfg.location = data.cfg.dip.pos;
ft_sourceplot(cfg, mri)

%%
% compute an averaged ERP over trials

cfg = [];
timelock = ft_timelockanalysis(cfg, data);

%%
% plot the averaged ERP

cfg = [];
cfg.layout = 'elec1010.lay';
ft_multiplotER(cfg, timelock);
