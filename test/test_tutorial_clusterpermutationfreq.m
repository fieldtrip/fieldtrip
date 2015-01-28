function test_tutorial_clusterpermutationfreq(dataset, datadir)

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_tutorial_eventrelatedstatistics
% TEST ft_freqanalysis ft_multiplotTFR ft_singleplotTFR ft_freqstatistics
% TEST ft_topoplotTFR ft_clustTFRplot ft_megplanar ft_combineplanar

global ft_default;
ft_default.feedback = 'no';

if nargin==0
  dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
  datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/cluster_permutation_freq');
end

%% PREprocessing
% find the interesting segments of data
cfg = [];                                    % empty configuration
cfg.dataset                 = dataset;       % name of CTF dataset  
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

% find the interesting segments of data
cfg = [];                                           % empty configuration
cfg.dataset                 = dataset;       % name of CTF dataset  
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 9;                    % event value of FC
cfg = ft_definetrial(cfg);            
  
% remove the trials that have artifacts from the trl
cfg.trl([2, 3, 4, 30, 39, 40, 41, 45, 46, 47, 51, 53, 59, 77, 85],:) = [];

% preprocess the data
cfg.channel   = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfg.demean    = 'yes';                              % do baseline correction with the complete trial

dataFC = ft_preprocessing(cfg);

%% planar gradient and TFA
cfg = [];
cfg.planarmethod = 'sincos';
% prepare_neighbours determines with what sensors the planar gradient is computed
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, dataFC);
  
dataFIC_planar = ft_megplanar(cfg, dataFIC);
dataFC_planar  = ft_megplanar(cfg, dataFC);

cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'MEG';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.foi        = 20;
cfg.toi        = [-1:0.05:2.0];
cfg.t_ftimwin  = 7./cfg.foi; %7 cycles
cfg.keeptrials = 'yes';

freqFC_planar  = ft_freqanalysis(cfg, dataFC_planar);
freqFIC_planar = ft_freqanalysis(cfg, dataFIC_planar);

cfg = [];
freqFIC_planar_cmb = ft_combineplanar(cfg, freqFIC_planar);
freqFC_planar_cmb = ft_combineplanar(cfg, freqFC_planar);
    
freqFIC_planar_cmb.grad = dataFIC.grad;
freqFC_planar_cmb.grad = dataFC.grad;

%% Permutation test
cfg = [];
cfg.channel          = {'MEG'};
cfg.latency          = 'all';
cfg.frequency        = 20;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, dataFC);

design = zeros(1,size(freqFIC_planar_cmb.powspctrm,1) + size(freqFC_planar_cmb.powspctrm,1));
design(1,1:size(freqFIC_planar_cmb.powspctrm,1)) = 1;
design(1,(size(freqFIC_planar_cmb.powspctrm,1)+1):(size(freqFIC_planar_cmb.powspctrm,1)+...
  size(freqFC_planar_cmb.powspctrm,1))) = 2;

cfg.design           = design;
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, freqFIC_planar_cmb, freqFC_planar_cmb);

%% Plotting the results
cfg = [];
cfg.jackknife = 'no';  
freqFIC_planar_cmb = ft_freqdescriptives(cfg, freqFIC_planar_cmb);
freqFC_planar_cmb  = ft_freqdescriptives(cfg, freqFC_planar_cmb);

stat.raweffect = freqFIC_planar_cmb.powspctrm - freqFC_planar_cmb.powspctrm;
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'raweffect';
cfg.zlim   = [-1e-27 1e-27];
cfg.layout = 'CTF151.lay';
ft_clusterplot(cfg, stat);


%% Prepropcessing and freqanalysis on planar data
cfg = [];
cfg.toilim = [-1.0 0];
dataFIC_baseline = ft_redefinetrial(cfg, dataFIC);

cfg = [];
cfg.toilim = [0.6 1.6];
dataFIC_activation = ft_redefinetrial(cfg, dataFIC);

dataFIC_baseline.time = dataFIC_activation.time;

cfg = [];
cfg.planarmethod = 'sincos';
% prepare_neighbours determines with what sensors the planar gradient is computed
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, dataFIC_activation);
dataFIC_activation_planar = ft_megplanar(cfg, dataFIC_activation);
dataFIC_baseline_planar   = ft_megplanar(cfg, dataFIC_baseline); 

cfg = [];
cfg.output = 'pow';
cfg.channel = 'MEG';
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 20;
cfg.toi = [0.6:0.05:1.6];
cfg.t_ftimwin = 7./cfg.foi; %7 cycles
cfg.keeptrials = 'yes';

freqFIC_activation_planar = ft_freqanalysis(cfg, dataFIC_activation_planar);
freqFIC_baseline_planar   = ft_freqanalysis(cfg, dataFIC_baseline_planar); 

cfg = [];
freqFIC_baseline_planar_cmb   = ft_combineplanar(cfg,freqFIC_baseline_planar);
freqFIC_activation_planar_cmb = ft_combineplanar(cfg,freqFIC_activation_planar);

freqFIC_baseline_planar_cmb.grad   = dataFIC.grad;
freqFIC_activation_planar_cmb.grad = dataFIC.grad;  

%% Permutation test
cfg = [];
cfg.channel          = {'MEG', '-MLP31', '-MLO12'};
cfg.channel          = {'MEG'};
cfg.latency          = [0.8 1.4];
cfg.method           = 'montecarlo';
cfg.frequency        = 20;
cfg.statistic        = 'ft_statfun_actvsblT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freqFIC_activation_planar_cmb);

ntrials = size(freqFIC_activation_planar_cmb.powspctrm,1);
design  = zeros(2,2*ntrials);
design(1,1:ntrials) = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials) = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];

cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;

[stat] = ft_freqstatistics(cfg, freqFIC_activation_planar_cmb, freqFIC_baseline_planar_cmb);

%% Plotting the results
cfg = [];
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.zlim   = [-4 4];
cfg.layout = 'CTF151.lay';
ft_clusterplot(cfg, stat);

%% Within subjects experiment

load(fullfile(datadir, 'GA_TFR_orig.mat'));

cfg = [];
cfg.channel          = {'MEG'};
cfg.latency          = [0 1.8];
cfg.frequency        = 20;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% specifies with which sensors other sensors can form clusters
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, GA_TFRFC);

subj = 10;
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;


failed = true;
for i=1:100 % this can fail sometime, like every second iteration or so, 100 is a *really* conservative number here
  [stat] = ft_freqstatistics(cfg, GA_TFRFIC, GA_TFRFC);
  
  % Plotting the results
  pcfg = [];
  pcfg.alpha  = 0.025;
  pcfg.parameter = 'stat';
  pcfg.zlim   = [-4 4];
  pcfg.layout = 'CTF151.lay';
  try
    ft_clusterplot(pcfg, stat);
    failed = false;
    break;
  end
end

if failed
  error('ft_clusterplot fails cause p was too high');
end
