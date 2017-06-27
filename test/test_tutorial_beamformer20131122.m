function test_tutorial_beamformer20131122

% MEM 10gb
% WALLTIME 02:30:00

% TEST ft_sourceanalysis ft_prepare_leadfield

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];
ft_default.feedback = 'no';

% this test script represents the MATLAB code from http://www.fieldtriptoolbox.org/tutorial/beamformer
% as downloaded from the wiki on 22 November 2013

cd(dccnpath('/home/common/matlab/fieldtrip/data'));

% find the interesting segments of data
cfg = [];                                           % empty configuration
cfg.dataset                 = 'Subject01.ds';       % name of CTF dataset
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 3;                    % event value of FIC
cfg = ft_definetrial(cfg);

cfg.trl([15, 36, 39, 42, 43, 49, 50, 81, 82, 84],:) = [];

cfg.channel   = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfg.demean    = 'yes';                              % do baseline correction with the complete trial

dataFIC = ft_preprocessing(cfg);

cfg = [];
cfg.toilim = [-0.5 0];
dataPre = ft_redefinetrial(cfg, dataFIC);

cfg.toilim = [0.8 1.3];
dataPost = ft_redefinetrial(cfg, dataFIC);

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

mri = ft_read_mri('Subject01.mri');
cfg = [];
cfg.write      = 'no';
cfg.coordsys   = 'ctf';
[segmentedmri] = ft_volumesegment(cfg, mri);

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);

cfg                 = [];
cfg.grad            = freqPost.grad;
cfg.vol             = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG','-MLP31', '-MLO12'};
cfg.grid.resolution = 1;   % use a 3-D grid with a 1 cm resolution
cfg.grid.unit       = 'cm';
[grid] = ft_prepare_leadfield(cfg);

cfg              = [];
cfg.method       = 'dics';
cfg.frequency    = 18;
cfg.grid         = grid;
cfg.vol          = vol;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = 0;

sourcePost_nocon = ft_sourceanalysis(cfg, freqPost);

mri = ft_read_mri('Subject01.mri');
mri = ft_volumereslice([], mri);

cfg            = [];
cfg.downsample = 2;
cfg.parameter = 'avg.pow';
sourcePostInt_nocon  = ft_sourceinterpolate(cfg, sourcePost_nocon , mri);

cfg              = [];
cfg.method       = 'slice';
cfg.funparameter = 'avg.pow';
figure
ft_sourceplot(cfg,sourcePostInt_nocon);

sourceNAI = sourcePost_nocon;
sourceNAI.avg.pow = sourcePost_nocon.avg.pow ./ sourcePost_nocon.avg.noise;

cfg = [];
cfg.downsample = 2;
cfg.parameter = 'avg.pow';
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

sourceDiff = sourcePost_con;
sourceDiff.avg.pow = (sourcePost_con.avg.pow - sourcePre_con.avg.pow) ./ sourcePre_con.avg.pow;

mri = ft_read_mri('Subject01.mri');
mri = ft_volumereslice([], mri);

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
ft_sourceplot(cfg, sourceDiffInt);

cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2];
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, sourceDiffInt);

cfg = [];
cfg.coordsys      = 'ctf';
cfg.nonlinear     = 'no';
sourceDiffIntNorm = ft_volumenormalise(cfg, sourceDiffInt);

cfg = [];
cfg.method        = 'ortho';
cfg.interactive   = 'yes';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2];
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, sourceDiffIntNorm);

cfg = [];
cfg.coordsys      = 'ctf';
cfg.nonlinear     = 'no';
sourceDiffIntNorm = ft_volumenormalise(cfg, sourceDiffInt);

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
ft_sourceplot(cfg, sourceDiffIntNorm);
view ([90 0])
