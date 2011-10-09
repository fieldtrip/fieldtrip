function test_ft_movieplotTFR

% TEST test_ft_movieplotTFR
% TEST ft_movieplotTFR ft_movieplotER

% the frequency analysis is based on the tutorials

load /home/common/matlab/fieldtrip/data/ftp/tutorial/timefrequencyanalysis/dataFIC.mat

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

figure
cfg = [];
cfg.layout = 'CTF151.lay';
ft_movieplotER(cfg, timelockFIC);

figure
cfg = [];
cfg.layout = 'CTF151.lay';
ft_movieplotTFR(cfg, freqFIC);

