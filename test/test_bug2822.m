function failed_bug2822

% WALLTIME 00:45:00
% MEM 6gb
% DEPENDENCY

% this is to test the implementation of the frequency domain MNE reconstruction

%% Start-up
%path_ft = '/home/electromag/lucamb/fieldtrip-dev';
%addpath(path_ft)

%% Load data (dataFIC)
path_to_load = dccnpath('/home/common/matlab/fieldtrip');
%load(dccnpath([path_to_load, '/data/test/dataFIC.mat']))

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

dataFIC = ft_preprocessing(cfg);

%% Frequency analysis

cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = 'MEG';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 18;
freq = ft_freqanalysis(cfg, dataFIC);

%% Genuine Headmodel and leadfield
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/sourcemodel/Subject01_headmodel.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/sourcemodel/Subject01_sourcemodel_15684.mat'));

headmodel = ft_convert_units(headmodel, 'cm'); % same as grad
sourcemodel = ft_convert_units(sourcemodel, 'cm'); % same as grad

%% Plot head model

% load vol                                       % volume conduction model
figure;
hold on;
ft_plot_headmodel(headmodel, 'facecolor', 'none');alpha 0.5;
ft_plot_mesh(sourcemodel, 'edgecolor', 'none'); camlight
ft_plot_sens(dataFIC.grad, 'style', '*b');

%% Precompute the leadfields
cfg = [];
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel;
cfg.singleshell.batchsize = 2500;
leadfield = ft_prepare_leadfield(cfg, freq);
leadfield.inside = true(size(leadfield.pos,1),1);

%% Source analysis MNE
cfg             = [];
cfg.method      = 'mne';
cfg.frequency   = 18;
cfg.sourcemodel = leadfield;
cfg.lambda      = 0.001;
source_freq_mne = ft_sourceanalysis(cfg, freq);

%% Source analysis RV

% cfg = [];
% cfg.method = 'rv';
% cfg.frequency = 18;
% cfg.sourcemodel = sourcemodel;
% cfg.headmodel = headmodel;
% cfg.lambda = 0.001;
% source_freq_rv = ft_sourceanalysis(cfg, freq);

%% Plot

figure,
ft_plot_mesh(sourcemodel, 'vertexcolor', source_freq_mne.avg.pow);
