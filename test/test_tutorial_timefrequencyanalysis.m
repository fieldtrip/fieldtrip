function test_tutorial_timefrequencyanalysis(datadir)

% MEM 2000mb
% WALLTIME 00:20:00

% TEST ft_freqanalysis ft_preprocessing ft_multiplotTFR ft_singleplotTFR

if nargin==0
  datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/timefrequencyanalysis');
end

load(fullfile(datadir, 'dataFIC.mat'));

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -0.5:0.05:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
TFRhann = ft_freqanalysis(cfg, dataFIC);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.zlim         = [-3e-27 3e-27];
cfg.showlabels   = 'yes';
cfg.layout       = 'CTF151.lay';
figure
ft_multiplotTFR(cfg, TFRhann);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.maskstyle    = 'saturation';
cfg.zlim         = [-3e-27 3e-27];
cfg.channel      = 'MRC15';
figure
ft_singleplotTFR(cfg, TFRhann);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.maskstyle    = 'saturation';
cfg.xlim         = [0.9 1.3];
cfg.zlim         = [-1.5e-27 1.5e-27];
cfg.ylim         = [15 20];
cfg.showlabels   = 'markers';
figure
ft_topoplotTFR(cfg, TFRhann);


cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MRC15';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:1:30;
cfg.t_ftimwin    = 7./cfg.foi;  % 7 cycles per time window
cfg.toi          = -0.5:0.05:1.5;
TFRhann7 = ft_freqanalysis(cfg, dataFIC);

cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.maskstyle    = 'saturation';
cfg.zlim         = [-3e-27 3e-27];
cfg.channel      = 'MRC15';
cfg.interactive  = 'no';
figure
ft_singleplotTFR(cfg, TFRhann7);

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MRC15';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:1:30;
cfg.t_ftimwin    = 4./cfg.foi;
cfg.toi          = -0.5:0.05:1.5;
TFRhann4 = ft_freqanalysis(cfg, dataFIC);

cfg.t_ftimwin    = 5./cfg.foi;
TFRhann5 = ft_freqanalysis(cfg, dataFIC);

cfg.t_ftimwin    = 10./cfg.foi;
TFRhann10 = ft_freqanalysis(cfg, dataFIC);

cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'MEG';
cfg.method     = 'mtmconvol';
cfg.foi        = 1:2:30;
cfg.t_ftimwin  = 5./cfg.foi;
cfg.tapsmofrq  = 0.4 *cfg.foi;
cfg.toi        = -0.5:0.05:1.5;
cfg.pad        = 'maxperlen';
TFRmult = ft_freqanalysis(cfg, dataFIC);


cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.zlim         = [-3e-27 3e-27];
cfg.showlabels   = 'yes';
cfg.layout       = 'CTF151.lay';
figure
ft_multiplotTFR(cfg, TFRmult)


cfg = [];
cfg.channel    = 'MEG';
cfg.method     = 'wavelet';
cfg.width      = 7;
cfg.output     = 'pow';
cfg.foi        = 1:2:30;
cfg.toi        = -0.5:0.05:1.5;
TFRwave = ft_freqanalysis(cfg, dataFIC);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.zlim         = [-3e-25 3e-25];
cfg.showlabels   = 'yes';
cfg.layout       = 'CTF151.lay';
figure
ft_multiplotTFR(cfg, TFRwave)
