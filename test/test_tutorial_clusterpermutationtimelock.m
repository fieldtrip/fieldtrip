function test_tutorial_clusterpermutationtimelock(dataset, datadir)

% MEM 2gb
% WALLTIME 00:10:00

% TEST test_tutorial_eventrelatedstatistics
% TEST ft_timelockanalysis ft_multiplotER ft_singleplotER ft_timelockstatistics
% TEST ft_topoplotER ft_clusterplot ft_megplanar ft_combineplanar

global ft_default;
ft_default.feedback = 'no';

if nargin==0
  dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
  datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/cluster_permutation_timelock');
end

%% Reading the FIC data
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

%% assert
% data = dataFIC_LP;
% load(fullfile(datadir, 'common', 'matlab', 'fieldtrip', 'data', 'ftp', 'tutorial', 'cluster_permutation_timelock', 'dataFIC_LP.mat'))
% data = rmfield(data, 'cfg');
% dataFIC_LP = rmfield(dataFIC_LP, 'cfg');
% data = rmfield(data, 'grad');
% dataFIC_LP = rmfield(dataFIC_LP, 'grad');
% assert(isequal(data, dataFIC_LP))

%% Reading the FC data
% find the interesting segments of data
cfg = [];                                           % empty configuration
cfg.dataset                 = dataset;       % name of CTF dataset  
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 9;                    % trigger value for fully incongruent (FC)
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

%% assert
% data = dataFC_LP;
% load(fullfile(datadir, 'common', 'matlab', 'fieldtrip', 'data', 'ftp', 'tutorial', 'cluster_permutation_timelock', 'dataFC_LP.mat'))
% data = rmfield(data, 'cfg');
% dataFC_LP = rmfield(dataFC_LP, 'cfg');
% data = rmfield(data, 'hdr');
% dataFC_LP = rmfield(dataFC_LP, 'hdr');
% data = rmfield(data, 'grad');
% dataFC_LP = rmfield(dataFC_LP, 'grad');
%assert(isequal(data, dataFC_LP))
% warning('This fails?'); % TODO check why this fails!

%% timelock analysis
cfg = [];
cfg.keeptrials = 'yes';
timelockFIC = ft_timelockanalysis(cfg, dataFIC_LP);
timelockFC  = ft_timelockanalysis(cfg, dataFC_LP);

%% Permutation test

cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to evaluate 
                                 % the effect at the sample level
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the permutation distribution. 
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is required for a selected 
%  cfg.neighbours = neighbours;  % see below for an explanation
                                 % sample to be included in the clustering algorithm (default=0).
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 100;      % number of draws from the permutation distribution

design = zeros(1,size(timelockFIC.trial,1) + size(timelockFC.trial,1));
design(1,1:size(timelockFIC.trial,1)) = 1;
design(1,(size(timelockFIC.trial,1)+1):(size(timelockFIC.trial,1) + size(timelockFC.trial,1)))= 2;

cfg.design = design;             % design matrix
cfg.ivar  = 1;                   % number or list with indices, independent variable(s)

cfg.channel = {'MEG'};          % cell-array with selected channel labels
cfg.latency = [0 1];            % time interval over which the experimental 
                                % conditions must be compared (in seconds)
cfg_neighb.method = 'distance';         
neighbours = ...            % specifies with which sensors other sensors
ft_prepare_neighbours(...       % can form clusters
cfg_neighb, dataFC_LP);
cfg.neighbours = neighbours;

[stat] = ft_timelockstatistics(cfg, timelockFIC, timelockFC);

%% Plotting the results
cfg = [];
cfg.keeptrials = 'no';  % now keep only the subject-wise average, not the single trials
avgFIC = ft_timelockanalysis(cfg, dataFIC_LP);
avgFC  = ft_timelockanalysis(cfg, dataFC_LP);

% Copy the entire timelockanalysis structure to preserve all
% the information it holds in addition to ERP averages.
raweffectFICvsFC     = avgFIC;
% Then take the difference of the averages.
raweffectFICvsFC.avg = avgFIC.avg - avgFC.avg;  

% Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
pos_cluster_pvals = [stat.posclusters(:).prob];
% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
neg = ismember(stat.negclusterslabelmat, neg_signif_clust);

pos = stat.posclusterslabelmat == 1;	% or == 2, or 3, etc.
neg = stat.negclusterslabelmat == 1;

timestep = 0.05;		% timestep between time windows for each subplot (in seconds)
sampling_rate = dataFC_LP.fsample;	% Data has a temporal resolution of 300 Hz
sample_count = length(stat.time);
					% number of temporal samples in the statistics object
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples

for k = 1:20;
     subplot(4,5,k);
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
     cfg.zlim = [-2.5e-13 2.5e-13];
   % If a channel reaches this significance, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).
   
   % Next, check which channels are significant over the
   % entire time interval of interest.
     pos_int = all(pos(:, m(k):m(k+1)), 2);
     neg_int = all(neg(:, m(k):m(k+1)), 2);

     cfg.highlight = 'on';
   % Get the index of each significant channel
     cfg.highlightchannel = find(pos_int | neg_int);
     cfg.comment = 'xlim';   
     cfg.commentpos = 'title';   
     cfg.layout = 'CTF151.lay';
     ft_topoplotER(cfg, raweffectFICvsFC);   
end

%% Using planar gradient data
cfg = [];
cfg.planarmethod = 'sincos';
cfg.neighbours = neighbours;
timelockFIC_planar = ft_megplanar(cfg, timelockFIC);
timelockFC_planar  = ft_megplanar(cfg, timelockFC);
       
timelockFIC_planar_cmb = ft_combineplanar(cfg, timelockFIC_planar);
timelockFC_planar_cmb  = ft_combineplanar(cfg, timelockFC_planar);
  
timelockFIC_planar_cmb.grad = timelockFIC.grad;  % add the gradiometer structure
timelockFC_planar_cmb.grad  = timelockFC.grad;

cfg = [];
cfg.channel = {'MEG'};
cfg.latency = [0 1];
cfg.neighbours = neighbours;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 100;

design = zeros(1,size(timelockFIC_planar_cmb.trial,1) + size(timelockFC_planar_cmb.trial,1));
design(1,1:size(timelockFIC_planar_cmb.trial,1)) = 1;
design(1,(size(timelockFIC_planar_cmb.trial,1)+1):(size(timelockFIC_planar_cmb.trial,1) + size(timelockFC_planar_cmb.trial,1)))= 2;

cfg.design = design;
cfg.ivar = 1;

[stat] = ft_timelockstatistics(cfg, timelockFIC_planar_cmb, timelockFC_planar_cmb);

%% assert here
%save stat_ERF_planar_FICvsFC stat

%% Plotting again

cfg = [];
cfg.keeptrials = 'no';   % now only the average, not the single trials
avgFIC_planar = ft_timelockanalysis(cfg, timelockFIC_planar);
avgFC_planar  = ft_timelockanalysis(cfg, timelockFC_planar);
cfg = [];
avgFIC_planar_cmb = ft_combineplanar(cfg, avgFIC_planar);
avgFC_planar_cmb  = ft_combineplanar(cfg, avgFC_planar);

% make a copy of the structure to contain the raw effect and all other information
% and subtract the averages from each other
raweffectFICvsFC     = avgFIC_planar_cmb;
raweffectFICvsFC.avg = avgFIC_planar_cmb.avg - avgFC_planar_cmb.avg;  

figure;  
timestep = 0.05;		%(in seconds)
sampling_rate = dataFC_LP.fsample;
sample_count = length(stat.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples

pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% Remember to do the same for negative clusters if you want them!

for k = 1:20;
     subplot(4,5,k);   
     cfg = [];
     cfg.xlim =[j(k) j(k+1)];
     cfg.zlim = [-1.0e-13 1.0e-13];   
     pos_int = all(pos(:, m(k):m(k+1)), 2);
     cfg.highlight = 'on';
     cfg.highlightchannel = find(pos_int);
     cfg.comment = 'xlim';
     cfg.commentpos = 'title';
     cfg.layout = 'CTF151.lay';
     ft_topoplotER(cfg, raweffectFICvsFC);
end


%% WIthin-subjects experiments

load(fullfile(datadir, 'GA_ERF_orig.mat'))
cfg = [];
cfg.channel = {'MEG'};
cfg.latency = [0 1];

cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 500;

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

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

[stat] = ft_timelockstatistics(cfg, GA_FIC, GA_FC);

%% assert
% save stat_ERF_planar_FICvsFC_GA stat

%% Plotting
GA_FICvsFC = GA_FIC;
GA_FICvsFC.avg = GA_FIC.avg - GA_FC.avg;

figure;  
timestep = 0.05;      %(in seconds)
sampling_rate = dataFC_LP.fsample;
sample_count = length(stat.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples

pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

for k = 1:20;
     subplot(4,5,k);   
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   
     cfg.zlim = [-1.0e-13 1.0e-13];   
     pos_int = all(pos(:, m(k):m(k+1)), 2);
     cfg.highlight = 'on';
     cfg.highlightchannel = find(pos_int);       
     cfg.comment = 'xlim';   
     cfg.commentpos = 'title';   
     ft_topoplotER(cfg, GA_FICvsFC);
end  
