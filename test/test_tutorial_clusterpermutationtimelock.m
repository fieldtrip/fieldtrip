function test_tutorial_clusterpermutationtimelock(dataset, datadir)

% MEM 3gb
% WALLTIME 00:10:00
% DEPENDENCY ft_timelockanalysis ft_multiplotER ft_singleplotER ft_timelockstatistics ft_topoplotER ft_clusterplot ft_megplanar ft_combineplanar

if nargin<1 || isempty(dataset)
  dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
end
if nargin<2
  datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/cluster_permutation_timelock');
end


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

cfg = [];
cfg.trials = data_all.trialinfo == 9;
dataFC_LP = ft_redefinetrial(cfg, data_all);

%save dataFIC_LP dataFIC_LP
%save dataFC_LP dataFC_LP

cfg = [];
cfg.keeptrials = 'yes';
timelockFIC    = ft_timelockanalysis(cfg, dataFIC_LP);
timelockFC     = ft_timelockanalysis(cfg, dataFC_LP);

%% ## Permutation test

cfg                  = [];
cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'indepsamplesT'; % use the independent samples T-statistic as a measure to
                                   % evaluate the effect at the sample level
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       % alpha level of the sample-specific test statistic that
                                   % will be used for thresholding
cfg.clusterstatistic = 'maxsum';   % test statistic that will be evaluated under the
                                   % permutation distribution.
cfg.minnbchan        = 2;          % minimum number of neighborhood channels that is
                                   % required for a selected sample to be included
                                   % in the clustering algorithm (default=0).
% cfg.neighbours     = neighbours; % see below
cfg.tail             = 0;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;      % alpha level of the permutation test
cfg.numrandomization = 100;        % number of draws from the permutation distribution

n_fc  = size(timelockFC.trial, 1);
n_fic = size(timelockFIC.trial, 1);

cfg.design           = [ones(1,n_fic), ones(1,n_fc)*2]; % design matrix
cfg.ivar             = 1; % number or list with indices indicating the independent variable(s)

cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, dataFC_LP);

cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
                                 % which other sensors it can form clusters
cfg.channel       = {'MEG'};     % cell-array with selected channel labels
cfg.latency       = [0 1];       % time interval over which the experimental
                                 % conditions must be compared (in seconds)
[stat] = ft_timelockstatistics(cfg, timelockFIC, timelockFC);

% Save the output to disk:
%save stat_ERF_axial_FICvsFC stat;

%% ## Plotting of the results
%
% To plot the results of the permutation test, we use the plotting function **[ft_topoplotER](https://github.com/fieldtrip/fieldtrip/blob/release/ft_topoplotER.m)**. In doing so, we will plot a topography of the difference between the two experimental conditions (FIC and FC). On top of that, and for each timestep of interest, we will highlight the sensors which are members of significant clusters. First, however, we must calculate the difference between conditions using **[ft_math](https://github.com/fieldtrip/fieldtrip/blob/release/ft_math.m)**.
%
cfg    = [];
avgFIC = ft_timelockanalysis(cfg, dataFIC_LP);
avgFC  = ft_timelockanalysis(cfg, dataFC_LP);

% Then take the difference of the averages using ft_math
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectFICvsFC = ft_math(cfg, avgFIC, avgFC);

% We then construct a boolean matrix indicating whether a channel/time point belongs to a cluster that we deem interesting to inspect. This matrix has size [Number_of_MEG_channels x Number_of_time_samples], like stat.posclusterslabelmat. We'll make two such matrices: one for positive clusters (named pos), and one for negative (neg). All (channel,time)-pairs belonging to the large clusters whose probability of occurrence is sufficiently low in relation to the associated randomization distribution of clusterstats will be coded in the new boolean matrix as 1, and all those that don't will be coded as 0.
%
%
pos_cluster_pvals = [stat.posclusters(:).prob];

% Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
% cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
% respectively
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(stat.posclusterslabelmat, pos_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_clust         = find(neg_cluster_pvals < 0.025);
neg               = ismember(stat.negclusterslabelmat, neg_clust);

% Alternatively, we can manually select which clusters we want to plot. If we only want to see the extext of the first (i.e. most significant) positive and negative clusters, for instance, we can do so as follows:
%
pos = stat.posclusterslabelmat == 1; % or == 2, or 3, etc.
neg = stat.negclusterslabelmat == 1;

% To plot a sequence of twenty topographic plots equally spaced between 0 and 1 second, we define the vector j of time steps. These time intervals correspond to the samples m in stat and in the variables pos and neg. m and j must, therefore, have the same length.
%
% To be sure that your sample-based time windows align with your time windows in seconds, check the following:
%
timestep      = 0.05; % timestep between time windows for each subplot (in seconds)
sampling_rate = dataFC_LP.fsample; % Data has a temporal resolution of 300 Hz
sample_count  = length(stat.time);
% number of temporal samples in the statistics object
j = [0:timestep:1]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in M/EEG samples

% To plot the data use the following for-loop:
%
% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffectFICvsFC.label, stat.label);

for k = 1:20
   subplot(4,5,k);
   cfg = [];
   cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
   cfg.zlim = [-2.5e-13 2.5e-13];
   % If a channel is in a to-be-plotted cluster, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).

   % Next, check which channels are in the clusters over the
   % entire time interval of interest.
   pos_int = zeros(numel(raweffectFICvsFC.label),1);
   neg_int = zeros(numel(raweffectFICvsFC.label),1);
   pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

   cfg.highlight   = 'on';
   % Get the index of the to-be-highlighted channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.comment     = 'xlim';
   cfg.commentpos  = 'title';
   cfg.layout      = 'CTF151_helmet.mat';
   cfg.interactive = 'no';
   cfg.figure      = 'gca';
   ft_topoplotER(cfg, raweffectFICvsFC);
end

%% ## Using planar gradient data
%
% To perform the permutation test using synthetic planar gradient data, the data must first be converted using the functions **[ft_megplanar](https://github.com/fieldtrip/fieldtrip/blob/release/ft_megplanar.m)** and **[ft_combineplanar](https://github.com/fieldtrip/fieldtrip/blob/release/ft_combineplanar.m)**. These functions were described in the tutorial on event-related fields. After running these functions, the statistical analysis using **[ft_timelockstatistics ](https://github.com/fieldtrip/fieldtrip/blob/release/ft_timelockstatistics.m)** involves the same configuration options as for ordinary event-related averages (no synthetic planar gradients). There is only one additional step, which is needed to add the gradiometer structure to one of the planar gradient data sets.
%
cfg = [];
cfg.planarmethod   = 'sincos';
cfg.neighbours     = neighbours; % also here, neighbouring sensors needs to be defined
timelockFIC_planar = ft_megplanar(cfg, timelockFIC);
timelockFC_planar  = ft_megplanar(cfg, timelockFC);

timelockFIC_planar_cmb = ft_combineplanar(cfg, timelockFIC_planar);
timelockFC_planar_cmb  = ft_combineplanar(cfg, timelockFC_planar);

timelockFIC_planar_cmb.grad = timelockFIC.grad;  % add the gradiometer structure
timelockFC_planar_cmb.grad  = timelockFC.grad;

% Having calculated synthetic planar gradient data, one can use the same configuration parameters as used for the analysis of the original data.
%
cfg                  = [];
cfg.channel          = {'MEG'};
cfg.latency          = [0 1];
cfg.neighbours       = neighbours;
cfg.method           = 'montecarlo';
cfg.statistic        = 'indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 100;

n_fc  = size(timelockFC_planar_cmb.trial, 1);
n_fic = size(timelockFIC_planar_cmb.trial, 1);

cfg.design           = [ones(1,n_fic), ones(1,n_fc)*2]; % design matrix
cfg.ivar             = 1; % number or list with indices indicating the independent variable(s)

cfg.ivar = 1;

[stat] = ft_timelockstatistics(cfg, timelockFIC_planar_cmb, timelockFC_planar_cmb);

%save stat_ERF_planar_FICvsFC stat

% The output can also be obtained from [ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/cluster_permutation_timelock/stat_ERF_planar_FICvsFC.mat](ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/cluster_permutation_timelock/stat_ERF_planar_FICvsFC.mat). If you need to reload the statistics output, use:
%
%load stat_ERF_planar_FICvsFC

% We now calculate the raw effect in the average with planar gradient data using the following configuration:
%
cfg = [];
cfg.keeptrials = 'no';   % now only the average, not the single trials
avgFIC_planar  = ft_timelockanalysis(cfg, timelockFIC_planar);
avgFC_planar   = ft_timelockanalysis(cfg, timelockFC_planar);
cfg = [];
avgFIC_planar_cmb = ft_combineplanar(cfg, avgFIC_planar);
avgFC_planar_cmb  = ft_combineplanar(cfg, avgFC_planar);

% subtract avgFC from avgFIC
cfg = [];
cfg.operation    = 'subtract';
cfg.parameter    = 'avg';
raweffectFICvsFC = ft_math(cfg, avgFIC_planar_cmb, avgFC_planar_cmb);

% Using the following configuration for **[ft_topoplotER](https://github.com/fieldtrip/fieldtrip/blob/release/ft_topoplotER.m)** we can plot the raw effect and highlight the channels contributing to the largest cluster
%
figure;
timestep      = 0.05; %(in seconds)
sampling_rate = dataFC_LP.fsample;
sample_count  = length(stat.time);
j = [0:timestep:1]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in M/EEG samples

pos_cluster_pvals = [stat.posclusters(:).prob];
pos_clust         = find(pos_cluster_pvals < 0.025);
pos               = ismember(stat.posclusterslabelmat, pos_clust);

% Remember to do the same for negative clusters if you want them!

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffectFICvsFC.label, stat.label);

for k = 1:20
   subplot(4,5,k);
   cfg = [];
   cfg.xlim =[j(k) j(k+1)];
   cfg.zlim = [-1.0e-13 1.0e-13];
   pos_int = zeros(numel(raweffectFICvsFC.label),1);
   pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   cfg.highlight = 'on';
   cfg.highlightchannel = find(pos_int);
   cfg.comment = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout = 'CTF151_helmet.mat';
   cfg.figure = 'gca';
   ft_topoplotER(cfg, raweffectFICvsFC);
end

%% # Within-subjects experiments
load(fullfile(datadir, 'ERF_orig.mat'));

% The configuration looks as follows:
cfg         = [];
cfg.channel = {'MEG'};
cfg.latency = [0 1];

cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;

Nsubj  = 10;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

[stat] = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});

%%save stat_ERF_planar_FICvsFC_GA stat

% The output can also be obtained from [stat_ERF_planar_FICvsFC_GA.mat](ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/cluster_permutation_timelock/stat_ERF_planar_FICvsFC_GA.mat). If you need to reload the statistics output, use:
%load stat_ERF_planar_FICvsFC_GA

%% ## Plotting the results
% load individual subject data
load(fullfile(datadir, 'ERF_orig.mat'));

% calculate the grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_FIC        = ft_timelockgrandaverage(cfg, allsubjFIC{:});
GA_FC         = ft_timelockgrandaverage(cfg, allsubjFC{:});

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_FICvsFC    = ft_math(cfg, GA_FIC, GA_FC);

figure;
% define parameters for plotting
timestep      = 0.05; %(in seconds)
sampling_rate = dataFIC_LP.fsample;
sample_count  = length(stat.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples

% get relevant values
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(stat.posclusterslabelmat, pos_clust);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(GA_FICvsFC.label, stat.label);

% plot
for k = 1:20
   subplot(4,5,k);
   cfg            = [];
   cfg.xlim       = [j(k) j(k+1)];
   cfg.zlim       = [-5e-14 5e-14];
   pos_int        = zeros(numel(GA_FICvsFC.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   cfg.highlight  = 'on';
   cfg.highlightchannel = find(pos_int);
   cfg.comment    = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout     = 'CTF151_helmet.mat';
   cfg.figure     = 'gca';
   ft_topoplotER(cfg, GA_FICvsFC);
end
