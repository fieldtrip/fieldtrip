function test_ft_movieplotTFR

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_movieplotTFR ft_movieplotER

% the frequency analysis is based on the tutorials

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/timefrequencyanalysis/dataFIC.mat'));

cfg              = [];
timelockFIC      = ft_timelockanalysis(cfg, dataFIC);

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:4:26;                         % changed from original
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -0.5:0.2:1.5;                   % changed from original
freqFIC          = ft_freqanalysis(cfg, dataFIC);

% non interactive timelock movie
figure
cfg = [];
cfg.samperframe  = 30;
cfg.framespersec = 7;
cfg.movierpt     = 3;
cfg.layout       = 'CTF151.lay';
ft_movieplotER(cfg, timelockFIC);

% interactive timelock movie
figure
cfg = [];
cfg.interactive = 'yes';
cfg.layout      = 'CTF151.lay';
ft_movieplotER(cfg, timelockFIC);

% interactive TFR movie
figure
cfg = [];
cfg.layout = 'CTF151.lay';
ft_movieplotTFR(cfg, freqFIC);

% non interactive TFR movie along frequencies
figure
cfg = [];
cfg.interactive = 'no';
cfg.movietime   = 1;
cfg.movierpt    = 3;
cfg.layout      = 'CTF151.lay';
ft_movieplotTFR(cfg, freqFIC);

% non interactive TFR movie along frequencies
figure
cfg = [];
cfg.interactive = 'no';
cfg.moviefreq   = 2;
cfg.movierpt    = 3;
cfg.layout = 'CTF151.lay';
ft_movieplotTFR(cfg, freqFIC);

% ensure that all figures are updated before XUnit starts to close the figures
drawnow
close all


drawnow

