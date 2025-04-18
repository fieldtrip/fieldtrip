function failed_bug2822

% WALLTIME 00:45:00
% MEM 4gb
% DEPENDENCY

% this is to test the implementation of the frequency domain MNE reconstruction

%% Start-up

%% Load data (dataFIC)

% find the interesting segments of data
cfg = [];
cfg.dataset                 = dccnpath('/project/3031000.02/external/download/test/ctf/Subject01.ds');       % name of CTF dataset
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

%% Fake Headmodel and leadfield

load(dccnpath('/home/common/matlab/fieldtrip/template/headmodel/standard_singleshell.mat'))

% use 'icosahedron' private function to generate the mash
mesh = [];
[mesh.pnt, mesh.tri] = icosahedron642;
mesh.pnt = 5*mesh.pnt - repmat([ 0 3 -1.5],size(mesh.pnt,1),1);

%% Plot head model

% load vol                                       % volume conduction model
figure;
hold on;
ft_plot_headmodel(vol, 'facecolor', 'none');alpha 0.5;
ft_plot_mesh(mesh, 'edgecolor', 'none'); camlight
ft_plot_sens(dataFIC.grad, 'style', '*b');

%% Source analysis MNE

cfg = [];
cfg.method = 'mne';
cfg.rawtrial = 'yes';
cfg.frequency = 18;
cfg.sourcemodel = mesh;
cfg.headmodel = vol;
cfg.lambda = 0.001;
source_freq_mne = ft_sourceanalysis(cfg, freq);

%% Source analysis RV

cfg = [];
cfg.method = 'rv';
cfg.frequency = 18;
cfg.sourcemodel = mesh;
cfg.headmodel = vol;
cfg.lambda = 0.001;
source_freq_rv = ft_sourceanalysis(cfg, freq);

%% Plot

figure,
ft_plot_mesh(mesh, 'vertexcolor', source_freq_mne.avg.pow);
