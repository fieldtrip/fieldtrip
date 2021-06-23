function test_beamformer_whitening

% WALLTIME 03:00:00
% MEM 8gb
% DEPENDENCY ft_denoise_prewhiten ft_sourceanalysis

% this function tests the functionality of running a beamformer analysis on
% spatially whitened data. TODO: at the moment - just to be safe - the mags
% and grads are whitened separately, and later concatenated. This requires
% for now some wizardry with the sensor arrays (also to be safe, to be
% honest), but it should be possible to have this implemented in
% ft_appendsens. Currently that does not work

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/original/rikhenson');
subj    = 15;

%% read the data from all separate runs
for run = 1 % the original uses 1:6, but that takes too long and too much storage
  trialdef = fullfile(datadir, sprintf('Sub%02d', subj), 'MEEG', 'Trials', sprintf('run_%02d_trldef.txt', run));
  dataset  = fullfile(datadir, sprintf('Sub%02d', subj), 'MEEG', sprintf('run_%02d_sss.fif', run));

  [begsample, endsample, offset, trialtype] = textread(trialdef, '%d%d%d%s');

  trialcode = nan(size(trialtype));
  trialcode(strcmp(trialtype, 'Famous'))      = 1;
  trialcode(strcmp(trialtype, 'Unfamiliar'))  = 2;
  trialcode(strcmp(trialtype, 'Scrambled'))   = 3;

  % construct the trial definition matrix, usually done with FT_DEFINETRIAL
  trl = [begsample(:) endsample(:) offset(:) trialcode(:)];

  cfg         = [];
  cfg.dataset = dataset;
  cfg.trl     = trl;

  % MEG specific settings
  cfg.channel = 'MEG';
  cfg.demean  = 'yes';
  cfg.coilaccuracy = 1;
  data = ft_preprocessing(cfg);
end

%% perform spatial whitening of the data
% split for now into mags and grads
cfg         = [];
cfg.channel = 'megmag';
data_mag    = ft_preprocessing(cfg, data);
cfg.channel = 'meggrad';
data_grad   = ft_preprocessing(cfg, data);

% extract the baseline covariance matrices
cfg         = [];
cfg.latency = [-inf 0];
cfg.covariance = 'yes';
cfg.covariancewindow = [-inf 0];
baseline_mag  = ft_timelockanalysis(cfg, data_mag);
baseline_grad = ft_timelockanalysis(cfg, data_grad);

% determine the kappa for the the whitening
[u,s,v] = svd(baseline_mag.cov);
figure;plot(log10(diag(s)),'o');title('singular values mags');
[u,s,v] = svd(baseline_grad.cov);
figure;plot(log10(diag(s)),'o');title('singular values grads');
% -> kappa = 71 seems to be a good enough value

% run the actual spatial whitening
cfg             = [];
cfg.kappa       = 71;
data_mag_white  = ft_denoise_prewhiten(cfg, data_mag,  baseline_mag);
data_grad_white = ft_denoise_prewhiten(cfg, data_grad, baseline_grad); 
data_white      = ft_appenddata([], data_mag_white, data_grad_white);

% combine the sens, this is a little bit of a hack, and needs to be
% addressed in ft_appendsens, if the coilpos is the same across sensor
% arrays, then the tra matrices can be appended
assert(isequal(data_grad_white.grad.coilpos,data_mag_white.grad.coilpos));
sens1 = data_mag_white.grad;
sens2 = data_grad_white.grad;
sens  = sens1;
fnames = {'chanori' 'chanpos' 'chantype' 'chanunit' 'label' 'tra'};
for k = 1:numel(fnames)
  sens.(fnames{k}) = cat(1,sens1.(fnames{k}), sens2.(fnames{k}));
end
fnames = {'labelnew' 'chantypenew' 'chanunitnew' 'tra'};
for k = 1:numel(fnames)
  sens.balance.prewhiten.(fnames{k}) = cat(1,sens1.balance.prewhiten.(fnames{k}), sens2.balance.prewhiten.(fnames{k}));
end
data_white.grad = sens;

%% create the geometric objects needed for the forward and inverse models
% load the original MRI
mri_orig = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/original/rikhenson/Sub15/T1/mprage.nii'));

% load the positions of the anatomical fiducials (as provided by Rik)
load(dccnpath('/home/common/matlab/fieldtrip/data/test/original/rikhenson/Sub15/T1/mri_fids.mat'));

% the location of fiducials is expressed in original MRI coordinates
% ft_volumerealign needs them in voxel coordinates
vox_fids = ft_warp_apply(inv(mri_orig.transform), mri_fids);

cfg              = [];
cfg.fiducial.nas = vox_fids(1,:);
cfg.fiducial.lpa = vox_fids(2,:);
cfg.fiducial.rpa = vox_fids(3,:);
cfg.coordsys  = 'neuromag';
mri_realigned = ft_volumerealign(cfg, mri_orig);

cfg        = [];
cfg.output = 'brain';
seg        = ft_volumesegment(cfg, mri_realigned);

cfg             = [];
cfg.method      = 'projectmesh';
cfg.numvertices = 2000;
cfg.tissue      = 'brain';
brain           = ft_prepare_mesh(cfg, seg);

%% make the volume conduction model
cfg        = [];
cfg.method = 'singleshell';
headmodel  = ft_prepare_headmodel(cfg, brain);

cfg = [];
cfg.resolution = 7;
cfg.headmodel = headmodel;  % from FT
%cfg.grad      = sens; % from FT
cfg.senstype  = 'meg';
cfg.singleshell.batchsize = 2000;
sourcemodel = ft_prepare_leadfield(cfg, data_white);

%% perform lcmv source reconstruction
cfg        = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 0];
cfg.covariance = 'yes';
tlck_white = ft_timelockanalysis(cfg, data_white);

[u,s,v] = svd(tlck_white.cov, 'econ');

cfg             = [];
cfg.headmodel   = headmodel;
cfg.grad        = sens;
cfg.sourcemodel = sourcemodel;
cfg.method      = 'lcmv';
cfg.lcmv.kappa        = 71;
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.fixedori     = 'yes';
cfg.lcmv.weightnorm   = 'unitnoisegain';
source = ft_sourceanalysis(cfg, tlck_white);

