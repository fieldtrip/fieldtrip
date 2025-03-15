function test_issue1473

% WALLTIME 00:45:00
% MEM 6gb
% DEPENDENCY

% this is to test the implementation of the MNE reconstruction, w.r.t.
% prewhitening

%% create a data structure 
% find the interesting segments of data
cfg = [];
cfg.dataset                 = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');       % name of CTF dataset
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 3;                    % event value of FIC
cfg = ft_definetrial(cfg);

% remove the trials that have artifacts from the trl
cfg.trl([15, 36, 39, 42, 43, 49, 50, 81, 82, 84],:) = [];

% preprocess the data
cfg.channel   = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfg.demean    = 'yes';                              % do baseline correction with the complete trial
cfg.baselinewindow = [-0.2 0];
cfg.bpfilter  = 'yes';
cfg.bpfreq    = [0.5 40];
cfg.bpfilttype = 'firws';
cfg.usefftfilt = 'yes';

data = ft_preprocessing(cfg);

%% prewhiten the data
cfg         = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-inf 0];
tlck        = ft_timelockanalysis(cfg, data);

data_white = ft_denoise_prewhiten([], data, tlck);
tlck_white = ft_timelockanalysis(cfg, data_white);


%% Genuine headmodel and sourcemodel
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/sourcemodel/Subject01_headmodel.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/sourcemodel/Subject01_sourcemodel_15684.mat'));

headmodel = ft_convert_units(headmodel, 'cm'); % same as grad
sourcemodel = ft_convert_units(sourcemodel, 'cm'); % same as grad
sourcemodel.inside = true(size(sourcemodel.pos,1),1);

%% Plot headmodel etc
figure;
hold on;
ft_plot_headmodel(headmodel, 'facecolor', 'none');alpha 0.5;
ft_plot_mesh(sourcemodel, 'edgecolor', 'none'); camlight
ft_plot_sens(data.grad, 'style', '*b');

%% Precompute the leadfields
cfg = [];
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel;
cfg.singleshell.batchsize = 2500;
leadfield        = ft_prepare_leadfield(cfg, data);

leadfield_white = ft_prepare_leadfield(cfg, data_white);



%% Source analysis MNE -> these give the same results (on september 11, 2020, using the master branch of FT)
cfg             = [];
cfg.method      = 'mne';
cfg.sourcemodel = leadfield;
cfg.mne.lambda  = 2;
cfg.mne.prewhiten = 'yes';
cfg.mne.scalesourcecov = 'yes';
source          = ft_sourceanalysis(cfg, tlck);

cfg.sourcemodel = leadfield_white;
source_white    = ft_sourceanalysis(cfg, tlck_white);

cfg.mne.prewhiten = 'no';
cfg.mne.lambda    = sqrt(cfg.mne.lambda); % the current state of affairs is an inconsistent treatment of lambda in the whitened, as compared to the non-whitened case, see issue 1307
source_white2     = ft_sourceanalysis(cfg, tlck_white);


