function test_tutorial_beamformer(datadir)

% MEM 8gb
% WALLTIME 03:30:00

% TEST ft_redefinetrial ft_freqanalysis ft_volumesegment ft_prepare_singleshell ft_sourceanalysis ft_prepare_leadfield ft_sourceinterpolate ft_sourceplot ft_volumenormalise

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];
ft_default.feedback = 'no';

if nargin==0
  % this is where the data should be located
  datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
end

load(fullfile(datadir, 'dataFIC.mat'));
load(fullfile(datadir, 'segmentedmri.mat'));
mri = ft_read_mri(fullfile(datadir, 'Subject01.mri'));

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

%try
  %if ~exist(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat'), 'file')
  cfg = [];
  cfg.write        = 'no';
  [segmentedmri] = ft_volumesegment(cfg, mri);
%catch
%  mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));
%  load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat'));
%end

%% Prepare head model
cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);

%% Prepare leadfield
cfg                 = [];
cfg.grad            = freqPost.grad;
cfg.vol             = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG','-MLP31', '-MLO12'};
cfg.grid.resolution = 1;   % use a 3-D grid with a 1 cm resolution
cfg.grid.unit       = 'cm';
[grid] = ft_prepare_leadfield(cfg);

%% Source analysis without contrasting condition

cfg              = []; 
cfg.method       = 'dics';
cfg.frequency    = 18;  
cfg.grid         = grid; 
cfg.vol          = vol;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = 0;

sourcePost_nocon = ft_sourceanalysis(cfg, freqPost);

%save sourcePost_nocon.mat sourcePost_nocon;

% call ft_volumereslice so that figures appear correct side up
mri = ft_volumereslice([], mri);

% Plot the result
cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'avg.pow';
sourcePostInt_nocon  = ft_sourceinterpolate(cfg, sourcePost_nocon , mri);

cfg              = [];
cfg.method       = 'slice';
cfg.funparameter = 'avg.pow';
figure
ft_sourceplot(cfg,sourcePostInt_nocon);

%% Compute and plot Neural Activity Index
sourceNAI = sourcePost_nocon;
sourceNAI.avg.pow = sourcePost_nocon.avg.pow ./ sourcePost_nocon.avg.noise;

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

%% Exercise 4: lead field normalization
cfg                 = [];
cfg.grad            = freqPost.grad;
cfg.vol             = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG','-MLP31', '-MLO12'};
cfg.grid.resolution = 1;   % use a 3-D grid with a 1 cm resolution
cfg.grid.unit       = 'cm';
cfg.normalize       = 'yes';
[gridn] = ft_prepare_leadfield(cfg);

cfg              = []; 
cfg.method       = 'dics';
cfg.frequency    = 18;  
cfg.grid         = gridn; 
cfg.vol          = vol;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = 0;
sourcePostn = ft_sourceanalysis(cfg, freqPost);

cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'avg.pow';
sourcePostIntn  = ft_sourceinterpolate(cfg, sourcePostn , mri);
cfg              = [];
cfg.method       = 'slice';
cfg.funparameter = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [.6 1.2]*1e-26;
cfg.opacitylim    = [.6 1.2]*1e-26;
cfg.opacitymap    = 'rampup';
figure
ft_sourceplot(cfg,sourcePostIntn);



%% Source analysis with constrasting condition
dataAll = ft_appenddata([], dataPre, dataPost);

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 4;
cfg.foilim    = [18 18];
freqAll = ft_freqanalysis(cfg, dataAll);

cfg              = [];
cfg.method       = 'dics';
cfg.frequency    = 18;
cfg.grid         = grid;
cfg.vol          = vol;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.realfilter   = 'yes';
sourceAll = ft_sourceanalysis(cfg, freqAll);

cfg.grid.filter = sourceAll.avg.filter;
sourcePre_con  = ft_sourceanalysis(cfg, freqPre );
sourcePost_con = ft_sourceanalysis(cfg, freqPost);

%save sourcePre_con sourcePre_con 
%save sourcePost_con sourcePost_con


% %% Plot results
% cfg            = [];
% cfg.downsample = 2;
% cfg.parameter  = 'avg.pow';
% sourcePostInt  = ft_sourceinterpolate(cfg, sourcePost , mri);
% 
% cfg              = [];
% cfg.method       = 'slice';
% cfg.funparameter = 'avg.pow';
% figure
% ft_sourceplot(cfg,sourcePostInt);
% 
% %% Plot Condition contrast

sourceDiff = sourcePost_con;
sourceDiff.avg.pow = (sourcePost_con.avg.pow - sourcePre_con.avg.pow) ./ sourcePre_con.avg.pow;

cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'avg.pow';
sourceDiffInt  = ft_sourceinterpolate(cfg, sourceDiff , mri);


cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0 1.2];
cfg.opacitylim    = [0 1.2];
cfg.opacitymap    = 'rampup';
figure
ft_sourceplot(cfg, sourceDiffInt);

%% Exercise 6: regularization
cfg              = [];
cfg.method       = 'dics';
cfg.frequency    = 18;
cfg.grid         = grid;
cfg.vol          = vol;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '0%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.realfilter   = 'yes';
sourceAll = ft_sourceanalysis(cfg, freqAll);
cfg.grid.filter = sourceAll.avg.filter;
source0Pre  = ft_sourceanalysis(cfg, freqPre );
source0Post = ft_sourceanalysis(cfg, freqPost);
cfg              = [];
cfg.method       = 'dics';
cfg.frequency    = 18;
cfg.grid         = grid;
cfg.vol          = vol;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '10%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.realfilter   = 'yes';
sourceAll = ft_sourceanalysis(cfg, freqAll);
cfg.grid.filter = sourceAll.avg.filter;
source10Pre  = ft_sourceanalysis(cfg, freqPre );
source10Post = ft_sourceanalysis(cfg, freqPost);

source0Diff = source0Post;
source0Diff.avg.pow = (source0Post.avg.pow - source0Pre.avg.pow) ./ source0Pre.avg.pow;
cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'avg.pow';
source0DiffInt  = ft_sourceinterpolate(cfg, source0Diff , mri);
cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0 1.2];
cfg.opacitylim    = [0 1.2];
cfg.opacitymap    = 'rampup';
figure
ft_sourceplot(cfg, source0DiffInt);

source10Diff = source10Post;
source10Diff.avg.pow = (source10Post.avg.pow - source10Pre.avg.pow) ./ source10Pre.avg.pow;
cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'avg.pow';
source10DiffInt  = ft_sourceinterpolate(cfg, source10Diff , mri);
cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0 1.2];
cfg.opacitylim    = [0 1.2];
cfg.opacitymap    = 'rampup';
figure
ft_sourceplot(cfg, source10DiffInt);


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
% cfg = [];
% cfg.nonlinear     = 'no';
% sourceDiffIntNorm = ft_volumenormalise(cfg, sourceDiffInt);


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
ft_sourceplot(cfg, sourceDiffIntNorm);
view ([90 0])


