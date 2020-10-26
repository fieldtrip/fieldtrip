function test_tutorial_beamformer_lcmv

% MEM 8gb
% WALLTIME 03:30:00
% DEPENDENCY ft_redefinetrial ft_freqanalysis ft_volumesegment ft_prepare_singleshell ft_sourceanalysis ft_prepare_leadfield ft_sourceinterpolate ft_sourceplot ft_volumenormalise

dataset = dccnpath('/home/common/matlab/fieldtrip/data/SubjectSEF.ds');

%% find the interesting segments of data
cfg                         = [];
cfg.dataset                 = dataset;
cfg.trialdef.eventtype      = 'Blue';
cfg.trialdef.prestim        = .1;        % .1 sec prior to trigger
cfg.trialdef.poststim       = .2;        % .2 sec following trigger
cfg.continuous              = 'yes';
cfg = ft_definetrial(cfg);

%% preprocess the data
cfg.channel                 = {'MEG'};
cfg.demean                  = 'yes';     % apply baselinecorrection
cfg.baselinewindow          = [-0.05 0]; % on basis of mean signal between -0.05 and 0
cfg.lpfilter                = 'yes';     % apply lowpass filter
cfg.lpfreq                  = 55;        % lowpass at 55 Hz
data = ft_preprocessing(cfg);

%% compute the covariance matrix
cfg                  = [];
cfg.covariance       = 'yes';
timelock             = ft_timelockanalysis(cfg, data);

figure; plot(timelock.time, timelock.avg)

%% view the results
cfg        = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.xlim   = [0.045 0.050];
ft_topoplotER(cfg, timelock);

%% calculate planar gradients
cfg                 = [];
cfg.feedback        = 'yes';
cfg.method          = 'template';
cfg.template        = 'ctf275_neighb.mat';
cfg.neighbours      = ft_prepare_neighbours(cfg, timelock);

cfg.planarmethod    = 'sincos';
timelock_planar     = ft_megplanar(cfg, timelock);

%% combine planar gradients
cfg                 = [];
timelock_planarcomb = ft_combineplanar(cfg, timelock_planar);


%% view the results
cfg        = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.xlim   = [0.045 0.050];
ft_topoplotER(cfg, timelock_planarcomb);

%% read and segment the subject's anatomical scan
load(dccnpath('/home/common/matlab/fieldtrip/data/SubjectSEF_mri.mat')); % matfile containing the realigned anatomical scan

cfg        = [];
cfg.output = 'brain';
seg = ft_volumesegment(cfg, mri);

%% make a figure of the mri and segmented volumes
segmentedmri           = seg;
segmentedmri.transform = mri.transform;
segmentedmri.anatomy   = mri.anatomy;
cfg                    = [];
cfg.funparameter       = 'brain';
ft_sourceplot(cfg, segmentedmri);

%% compute the subject's headmodel/volume conductor model
cfg                = [];
cfg.method         = 'singleshell';
headmodel          = ft_prepare_headmodel(cfg, seg);
headmodel          = ft_convert_units(headmodel, 'cm'); % mm to cm

%% create the subject specific grid
grad = ft_read_sens(dataset, 'senstype', 'meg');

cfg             = [];
cfg.grad        = grad;
cfg.headmodel   = headmodel;
cfg.resolution  = 0.5;
cfg.inwardshift = -1;
sourcemodel     = ft_prepare_sourcemodel(cfg);

%% make a figure of the single subject headmodel, and grid positions
figure;
ft_plot_sens(grad, 'style', '*b');
ft_plot_headmodel(headmodel, 'edgecolor', 'none'); alpha 0.4;
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));

%% create leadfield
cfg                  = [];
cfg.grad             = grad;  % gradiometer distances
cfg.headmodel        = headmodel;   % volume conduction headmodel
cfg.sourcemodel      = sourcemodel;
cfg.channel          = {'MEG'};
cfg.singleshell.batchsize = 2000;
lf                   = ft_prepare_leadfield(cfg);

%% create spatial filter using the lcmv beamformer
cfg                  = [];
cfg.method           = 'lcmv';
cfg.sourcemodel      = lf; % leadfield
cfg.headmodel        = headmodel; % volume conduction model (headmodel)
cfg.lcmv.keepfilter  = 'yes';
cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
source               = ft_sourceanalysis(cfg, timelock);

