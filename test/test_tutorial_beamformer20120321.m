function test_tutorial_beamformer20120321

% MEM 10gb
% WALLTIME 02:30:00

% TEST test_tutorial_beamformer
% TEST ft_redefinetrial ft_freqanalysis ft_volumesegment ft_prepare_singleshell ft_sourceanalysis ft_prepare_leadfield ft_sourceinterpolate ft_sourceplot ft_volumenormalise

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];
ft_default.feedback = 'no';

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/dataFIC.mat'));

%% Preprocess time windows of interest

cfg = [];
cfg.toilim = [-0.5 0];
dataPre = ft_redefinetrial(cfg, dataFIC);

cfg.toilim = [0.8 1.3];
dataPost = ft_redefinetrial(cfg, dataFIC);

%% Cross-spectral density

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 4;
cfg.foilim    = [18 18];
freqPre = ft_freqanalysis(cfg, dataPre);

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 4;
cfg.foilim    = [18 18];
freqPost = ft_freqanalysis(cfg, dataPost);

%% Compute (or load) the forward model)

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));
if ~exist(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat'), 'file')
  % segment the anatomical MRI
  cfg = [];
  cfg.write        = 'no';
  cfg.coordsys     = 'ctf';
  [segmentedmri] = ft_volumesegment(cfg, mri);
else
  % use the segmented MRI that is available
  load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat'));
end

%% Prepare head model
cfg = [];
vol = ft_prepare_singleshell(cfg, segmentedmri);

%% Prepare leadfield
cfg                 = [];
cfg.grad            = freqPre.grad;
cfg.vol             = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG','-MLP31', '-MLO12'};
cfg.grid.resolution = 1;   % use a 3-D grid with a 1 cm resolution
cfg.grid.unit = 'cm';
[grid] = ft_prepare_leadfield(cfg);

%% Source analysis
cfg              = [];
cfg.frequency    = 18;
cfg.method       = 'dics';
cfg.projectnoise = 'yes';
cfg.grid         = grid;
cfg.vol          = vol;
cfg.lambda       = 0;

sourcePre  = ft_sourceanalysis(cfg, freqPre );
sourcePost = ft_sourceanalysis(cfg, freqPost);

%save source sourcePre sourcePost


%% Plot results
cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'avg.pow';
sourcePostInt  = ft_sourceinterpolate(cfg, sourcePost , mri);

cfg              = [];
cfg.method       = 'slice';
cfg.funparameter = 'avg.pow';
figure
ft_sourceplot(cfg,sourcePostInt);

%% Plot Condition contrast

sourceDiff = sourcePost;
sourceDiff.avg.pow = (sourcePost.avg.pow - sourcePre.avg.pow) ./ sourcePre.avg.pow;

cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'avg.pow';
sourceDiffInt  = ft_sourceinterpolate(cfg, sourceDiff , mri);


cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2];
cfg.opacitymap    = 'rampup';
figure
ft_sourceplot(cfg, sourceDiffInt);


%% Orthogonal cut

cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2];
cfg.opacitymap    = 'rampup';
figure
ft_sourceplot(cfg, sourceDiffInt);


%% Template MRI

cfg = [];
cfg.coordsys      = 'ctf';
cfg.nonlinear     = 'no';
sourceDiffIntNorm = ft_volumenormalise(cfg, sourceDiffInt);

cfg = [];
cfg.method        = 'ortho';
%cfg.interactive   = 'yes';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2];
cfg.opacitymap    = 'rampup';
figure
ft_sourceplot(cfg, sourceDiffIntNorm);

%% Project to a surface
cfg = [];
cfg.coordsys   = 'ctf';
cfg.nonlinear  = 'no';
sourceDiffIntN = ft_volumenormalise(cfg, sourceDiffInt);


cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'avg.pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolorlim    = [0.0 1.2];
cfg.funcolormap    = 'jet';
cfg.opacitylim     = [0.0 1.2];
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'surface_white_both.mat';
cfg.surfdownsample = 10;
figure
ft_sourceplot(cfg, sourceDiffIntN);
view ([90 0])


%% Neural Activity Index (NAI)

sourceNAI = sourcePost;
sourceNAI.avg.pow = sourcePost.avg.pow ./ sourcePost.avg.noise;

cfg = [];
cfg.downsample = 2;
cfg.parameter  = 'avg.pow';
sourceNAIInt = ft_sourceinterpolate(cfg, sourceNAI , mri);


cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [4.0 6.2];
cfg.opacitylim    = [4.0 6.2];
cfg.opacitymap    = 'rampup';
figure
ft_sourceplot(cfg, sourceNAIInt);
