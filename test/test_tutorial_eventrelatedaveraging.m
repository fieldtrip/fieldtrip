function test_tutorial_eventrelatedaveraging(dataset)

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preprocessing ft_timelockanalysis ft_multiplotER ft_singleplotER ft_topoplotER ft_megplanar ft_combineplanar

% see http://www.fieldtriptoolbox.org/tutorial/eventrelatedaveraging
% this testscript corresponds to the version on the wiki at 23 December 2012

if nargin<1
  % use the tutorial dataset from home/common
  dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
else
  % use the dataset specified in the input, but check that it is called
  % Subject01.ds
  [p,n,e] = fileparts(dataset);
  if ~strcmp(n, 'Subject01')
    error('the dataset should be Subject01.ds');
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing

% find the interesting segments of data
cfg = [];                                           % empty configuration
cfg.dataset                 = dataset;       % name of CTF dataset
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 3;                    % trigger value for fully incongruent (FIC)
cfg = ft_definetrial(cfg);

% remove the trials that have artifacts from the trl
cfg.trl([15, 36, 39, 42, 43, 49, 50, 81, 82, 84],:) = [];

% preprocess the data
cfg.channel    = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfg.demean     = 'yes';
cfg.baselinewindow  = [-0.2 0];
cfg.lpfilter   = 'yes';                              % apply lowpass filter
cfg.lpfreq     = 35;                                 % lowpass at 35 Hz.

dataFIC_LP = ft_preprocessing(cfg);

plot(dataFIC_LP.time{1}, dataFIC_LP.trial{1}(130,:))

% find the interesting segments of data
cfg = [];                                           % empty configuration
cfg.dataset                 = dataset;       % name of CTF dataset
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 9;                    % trigger value for fully congruent (FC)
cfg = ft_definetrial(cfg);

% remove the trials that have artifacts from the trl
cfg.trl([2, 3, 4, 30, 39, 40, 41, 45, 46, 47, 51, 53, 59, 77, 85],:) = [];

% preprocess the data
cfg.channel    = {'MEG', '-MLP31', '-MLO12'};       % read all MEG channels except MLP31 and MLO12
cfg.demean     = 'yes';
cfg.baselinewindow  = [-0.2 0];
cfg.lpfilter   = 'yes';                              % apply lowpass filter
cfg.lpfreq     = 35;                                 % lowpass at 35 Hz.

dataFC_LP = ft_preprocessing(cfg);

% find the interesting segments of data
cfg = [];                                           % empty configuration
cfg.dataset                 = dataset;       % name of CTF dataset
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 5;                    % trigger value for initially congruent (IC)
cfg = ft_definetrial(cfg);

% remove the trials that have artifacts from the trl
cfg.trl([1, 2, 3, 4, 14, 15, 16, 17, 20, 35, 39, 40, 47, 78, 79, 80, 86],:) = [];

% preprocess the data
cfg.channel    = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfg.demean     = 'yes';
cfg.baselinewindow  = [-0.2 0];
cfg.lpfilter   = 'yes';                              % apply lowpass filter
cfg.lpfreq     = 35;                                 % lowpass at 35 Hz.

dataIC_LP = ft_preprocessing(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Timelockanalysis

cfg = [];
avgFIC = ft_timelockanalysis(cfg, dataFIC_LP);
avgFC = ft_timelockanalysis(cfg, dataFC_LP);
avgIC = ft_timelockanalysis(cfg, dataIC_LP);

cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 6;
cfg.layout = 'CTF151.lay';
cfg.ylim = [-3e-13 3e-13];
ft_multiplotER(cfg, avgFIC);

cfg = [];
cfg.showlabels = 'no';
cfg.fontsize = 6;
cfg.layout = 'CTF151.lay';
cfg.baseline = [-0.2 0];
cfg.xlim = [-0.2 1.0];
cfg.ylim = [-3e-13 3e-13];
ft_multiplotER(cfg, avgFC, avgIC, avgFIC);

cfg.xlim = [-0.2 1.0];
cfg.ylim = [-1e-13 3e-13];
cfg.channel = 'MLC24';
clf;
ft_singleplotER(cfg,avgFC, avgIC, avgFIC);

cfg = [];
cfg.xlim = [0.3 0.5];
cfg.colorbar = 'yes';
ft_topoplotER(cfg,avgFIC);

cfg = [];
cfg.xlim = [-0.2 : 0.1 : 1.0];  % Define 12 time intervals
cfg.zlim = [-2e-13 2e-13];      % Set the 'color' limits.
clf;
ft_topoplotER(cfg,avgFIC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the planar gradient

cfg                 = [];
cfg.feedback        = 'yes';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, avgFIC);

cfg.planarmethod    = 'sincos';
avgFICplanar        = ft_megplanar(cfg, avgFIC);

cfg = [];
avgFICplanarComb = ft_combineplanar(cfg,avgFICplanar);

cfg = [];
clf;
subplot(121);
cfg.xlim = [0.3 0.5];
cfg.zlim = 'maxmin';
cfg.colorbar = 'yes';
cfg.layout = 'CTF151.lay';
ft_topoplotER(cfg,avgFIC)
colorbar;
subplot(122);
cfg.zlim = 'maxabs';
cfg.layout = 'CTF151.lay';
ft_topoplotER(cfg,avgFICplanarComb);

