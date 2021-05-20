function test_pull1573(datadir)

% MEM 8gb
% WALLTIME 03:30:00
% DEPENDENCY ft_redefinetrial ft_freqanalysis ft_volumesegment ft_prepare_singleshell ft_sourceanalysis ft_prepare_leadfield ft_sourceinterpolate ft_sourceplot ft_volumenormalise

if nargin==0
  % datadir points by default to the fileserver at the Donders, use an
  % input argument if you run this function locally
  datadir = dccnpath('/home/common/matlab/fieldtrip/data/');
end

%% INITIALIZE
% this should be done outside the function, if it's to be run as a test
% function on a more or less daily basis

%% READ DATA AND PREPROCESS
% data has been got from;
% ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/SubjectSEF.zip
% put subject folder in working directory

cfg = [];
cfg.dataset = fullfile(datadir, 'SubjectSEF.ds');
cfg.trialdef.eventtype = 'Blue';
cfg.trialdef.prestim = 0.100;
cfg.trialdef.poststim = 0.200;
cfg.continuous = 'yes';

cfg = ft_definetrial(cfg);

cfg.channel = 'MEG';
cfg.demean  = 'yes';
cfg.baselinewindow = [-0.100 -0.010];

data = ft_preprocessing(cfg);

%% AVERAGING AND COMPUTATION OF THE COVARIANCE MATRIX

cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-inf 0];
baseline = ft_timelockanalysis(cfg, data);

cfg.covariancewindow = [-inf inf];
timelock = ft_timelockanalysis(cfg, data);

%% CREATE A WHITENED VERSION OF THE DATA
% for this the whitener should be different than the data that is to be
% used for the beamformer. Reason: one wants to whiten the noise, but keep
% the structure in the 'signal'. If you whiten with timelock, you will
% trivially get an Identity matrix as the data covariance.

cfg = [];
white_timelock = ft_denoise_prewhiten(cfg, timelock, baseline);

%% I guess the beamformer needs the original covariance?? NO, this is not what you needtest
%white_timelock.cov = timelock.cov; 

% this is also needed for interpolation later on
load(fullfile(datadir, 'SubjectSEF_mri.mat'));

%% try to skip the headmodel creation part, by loading the file from disk
try
  load(fullfile(datadir, 'ftp/tutorial/beamformer_lcmv', 'headmodel.mat'));
catch
  %% SEGMENT MRI
  cfg        = [];
  cfg.output = 'brain';
  segmented_mri = ft_volumesegment(cfg, mri);
  
  %% HEAD MODEL
  
  cfg        = [];
  cfg.method = 'singleshell';
  cfg.unit   = 'cm';
  headmodel  = ft_prepare_headmodel(cfg, segmented_mri);
end

%% SOURCE MODEL

cfg = [];
cfg.grad = timelock.grad;
cfg.headmodel = headmodel;
cfg.resolution = .75;
cfg.inwardshift = -.75;

sourcemodel = ft_prepare_sourcemodel(cfg);

%% CHECK COREGISTRATION

figure
ft_plot_sens(timelock.grad, 'style', '*b');
ft_plot_headmodel(headmodel, 'edgecolor', 'none');
alpha 0.4;
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside, :));

%% LEADFIELD

cfg = [];
cfg.grad = timelock.grad;
cfg.headmodel = headmodel;
cfg.sourcemodel = sourcemodel;
cfg.channel = 'MEG';
cfg.singleshell.batchsize = 1000;

leadfield = ft_prepare_leadfield(cfg);

% here's the trick, when using whitened data, the leadfield should also be
% based on the whitened data
cfg.grad = white_timelock.grad;
white_leadfield = ft_prepare_leadfield(cfg);

%% SOURCE ANALYSIS

n_components = 10; %% we can try different numbers

% "normal" lcmv
cfg = [];
cfg.method = 'lcmv';
cfg.sourcemodel = leadfield;
cfg.headmodel = headmodel;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.projectnoise = 'yes';
%cfg.lcmv.weightnorm = 'arraygain'; % this is not ideal for a one-to-one
%comparison, let's switch off for now.

cfg.lcmv.eigenspace = 'no';
source              = ft_sourceanalysis(cfg, timelock);

cfg.sourcemodel = white_leadfield;
source_white    = ft_sourceanalysis(cfg, white_timelock);

% eigenspace
cfg.sourcemodel         = leadfield;
cfg.lcmv.eigenspace     = n_components;
source_eigenspace       = ft_sourceanalysis(cfg, timelock);

cfg.sourcemodel         = white_leadfield;
cfg.lcmv.prewhitened    = 'yes';
source_white_eigenspace = ft_sourceanalysis(cfg, white_timelock);

% subspace
cfg.sourcemodel       = leadfield;
cfg.lcmv.eigenspace   = 'no';
cfg.lcmv.prewhitened  = 'no';
cfg.lcmv.subspace     = n_components;
source_subspace       = ft_sourceanalysis(cfg, timelock);

cfg.sourcemodel       = white_leadfield;
source_white_subspace = ft_sourceanalysis(cfg, white_timelock);

sources = {source source_white ...
           source_eigenspace source_white_eigenspace ...
           source_subspace source_white_subspace};
         
%% PLOT MOMENT AT 46 MS

close all

n_analyses = length(sources);

for analysis_index = 1:n_analyses
  
  source = sources{analysis_index};
  
  mom  = cat(1, source.avg.mom{:});
  bmom = std(mom(:,1:nearest(source.time,0)),[],2);
  sel  = nearest(source.time, [0.04 0.055]);
  pow  = abs(mean(mom(:,sel(1):sel(2)),2))./bmom;
  
  M{analysis_index} = abs(mom)./repmat(bmom, [1 numel(source.time)]);
  tmp = zeros(size(sourcemodel.inside));
  tmp(sourcemodel.inside) = pow;
  sourcemodel.pow = tmp;
%   
%   % pick time
%   cfg = [];
%   cfg.latency = 0.046;
%   
%   source = ft_selectdata(cfg, source);
%   
%   % I'm sure there's a nicer way to do this
%   mom_array = zeros(size(source.mom));
%   n_sources = length(mom_array);
%   
%   for source_index = 1:n_sources
%     moment = source.mom{source_index};
%     if isempty(moment)
%       mom_array(source_index) = NaN;
%     else
%       mom_array(source_index) = moment;
%     end
%   end
%   
%   mom_array = abs(mom_array);
%   
%   source.mom_array = mom_array;
  
  cfg = [];
  cfg.parameter = 'pow';
  %cfg.parameter = 'mom_array';
  
  source = ft_sourceinterpolate(cfg, sourcemodel, mri);
  %source = ft_sourceinterpolate(cfg, source, mri);
  
  cfg = [];
  cfg.funparameter = 'pow';
  cfg.funcolormap = 'viridis';
  %cfg.funparameter = 'mom_array';
  
  ft_sourceplot(cfg, source);
  
end

keyboard