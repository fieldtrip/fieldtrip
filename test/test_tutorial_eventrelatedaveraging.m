function test_tutorial_eventrelatedaveraging(dataset)

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preprocessing ft_timelockanalysis ft_multiplotER ft_singleplotER ft_topoplotER ft_megplanar ft_combineplanar

% see http://www.fieldtriptoolbox.org/tutorial/eventrelatedaveraging
% this testscript corresponds to the version on the wiki at 23 December 2012

if nargin<1
  % use the tutorial dataset from home/common
  dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
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

cfg                         = [];
cfg.dataset                 = dataset;
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.eventvalue     = [3 5 9]; % the values of the stimulus trigger for the three conditions
% 3 = fully incongruent (FIC), 5 = initially congruent (IC), 9 = fully congruent (FC)
cfg.trialdef.prestim        = 1; % in seconds
cfg.trialdef.poststim       = 2; % in seconds

cfg = ft_definetrial(cfg);

% remove the trials that have artifacts from the trl
cfg.trl([2, 5, 6, 8, 9, 10, 12, 39, 43, 46, 49, 52, 58, 84, 102, 107, 114, 115, 116, 119, 121, 123, 126, 127, 128, 133, 137, 143, 144, 147, 149, 158, 181, 229, 230, 233, 241, 243, 245, 250, 254, 260],:) = [];

% preprocess the data
cfg.channel    = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfg.demean     = 'yes';
cfg.baselinewindow  = [-0.2 0];
cfg.lpfilter   = 'yes';                              % apply lowpass filter
cfg.lpfreq     = 35;                                 % lowpass at 35 Hz.

data_all = ft_preprocessing(cfg);



cfg        = [];
cfg.trials = data_all.trialinfo == 3;
dataFIC_LP = ft_redefinetrial(cfg, data_all);

cfg        = [];
cfg.trials = data_all.trialinfo == 5;
dataIC_LP = ft_redefinetrial(cfg, data_all);

cfg = [];
cfg.trials = data_all.trialinfo == 9;
dataFC_LP = ft_redefinetrial(cfg, data_all);

%save dataFIC_LP dataFIC_LP
%save dataFC_LP dataFC_LP
%save dataIC_LP dataIC_LP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Timelockanalysis

cfg = [];
avgFIC = ft_timelockanalysis(cfg, dataFIC_LP);
avgFC = ft_timelockanalysis(cfg, dataFC_LP);
avgIC = ft_timelockanalysis(cfg, dataIC_LP);

cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 6;
cfg.layout = 'CTF151_helmet.mat';
cfg.ylim = [-3e-13 3e-13];
ft_multiplotER(cfg, avgFIC);

cfg = [];
cfg.showlabels = 'no';
cfg.fontsize = 6;
cfg.layout = 'CTF151_helmet.mat';
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
cfg.senstype        = 'meg';    
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
cfg.layout = 'CTF151_helmet.mat';
ft_topoplotER(cfg,avgFIC)
colorbar;
subplot(122);
cfg.zlim = 'maxabs';
cfg.layout = 'CTF151_helmet.mat';
ft_topoplotER(cfg,avgFICplanarComb);

