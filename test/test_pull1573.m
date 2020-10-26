%% INITIALIZE

clear variables
restoredefaultpath

% add fieldtrip <insert your path>
addpath /home/lau/matlab/fieldtrip/
ft_defaults

%% READ DATA AND PREPROCESS
% data has been got from;
% ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/SubjectSEF.zip
% put subject folder in working directory

cfg = [];
cfg.dataset = 'SubjectSEF.ds';
cfg.trialdef.eventtype = 'Blue';
cfg.trialdef.prestim = 0.100;
cfg.trialdef.poststim = 0.200;
cfg.continuous = 'yes';

cfg = ft_definetrial(cfg);

cfg.channel = 'MEG';
cfg.demean = 'yes';
cfg.baselinewindow = [-0.100 -0.010];

data = ft_preprocessing(cfg);

%% AVERAGING AND COMPUTATION OF THE COVARIANCE MATRIX

cfg = [];
cfg.covariance = 'yes';

timelock = ft_timelockanalysis(cfg, data);

%% CREATE A WHITENED VERSION

cfg = [];

white_timelock = ft_denoise_prewhiten(cfg, timelock, timelock);

% I guess the beamformer needs the original covariance??
white_timelock.cov = timelock.cov; 

%% SEGMENT MRI

load('SubjectSEF_mri.mat');

cfg = [];
cfg.output = 'brain';

segmented_mri = ft_volumesegment(cfg, mri);

%% HEAD MODEL

cfg = [];
cfg.method = 'singleshell';
cfg.unit = 'cm';

headmodel = ft_prepare_headmodel(cfg, segmented_mri);

%% SOURCE MODEL

cfg = [];
cfg.grad = timelock.grad;
cfg.headmodel = headmodel;
cfg.resolution = 1;
cfg.inwardshift = -1;

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

leadfield = ft_prepare_leadfield(cfg);

%% SOURCE ANALYSIS

n_components = 10; %% we can try different numbers

% "normal" lcmv
cfg = [];
cfg.method = 'lcmv';
cfg.sourcemodel = leadfield;
cfg.headmodel = headmodel;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.weightnorm = 'arraygain';
cfg.lcmv.eigenspace = 'no';

source = ft_sourceanalysis(cfg, timelock);
source_white = ft_sourceanalysis(cfg, white_timelock);

% eigenspace

cfg.lcmv.eigenspace = n_components;

source_eigenspace = ft_sourceanalysis(cfg, timelock);

cfg.lcmv.prewhitened = 'yes';

source_white_eigenspace = ft_sourceanalysis(cfg, white_timelock);

% subspace

cfg.lcmv.eigenspace = 'no';
cfg.lcmv.prewhitened = 'yes';
cfg.lcmv.subspace = n_components;

source_subspace = ft_sourceanalysis(cfg, timelock);
source_white_subspace = ft_sourceanalysis(cfg, white_timelock);

sources = {source source_white ...
           source_eigenspace source_white_eigenspace ...
           source_subspace source_white_subspace};
         
%% PLOT MOMENT AT 46 MS

close all

n_analyses = length(sources);

for analysis_index = 1:n_analyses
  
  source = sources{analysis_index};
  
  % pick time
  cfg = [];
  cfg.latency = 0.046;
  
  source = ft_selectdata(cfg, source);
  
  % I'm sure there's a nicer way to do this
  mom_array = zeros(size(source.mom));
  n_sources = length(mom_array);
  
  for source_index = 1:n_sources
    moment = source.mom{source_index};
    if isempty(moment)
      mom_array(source_index) = NaN;
    else
      mom_array(source_index) = moment;
    end
  end
  
  mom_array = abs(mom_array);
  
  source.mom_array = mom_array;
  
  cfg = [];
  cfg.parameter = 'mom_array';
  
  source = ft_sourceinterpolate(cfg, source, mri);
  
  cfg = [];
  cfg.funparameter = 'mom_array';
  
  ft_sourceplot(cfg, source);
  
end