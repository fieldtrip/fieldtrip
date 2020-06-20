function test_pull1456

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_sourceanalysis ft_eloreta

load(dccnpath('/home/common/matlab/fieldtrip/data/test/test_pull1456.mat'));

% Compute spectrum
cfg                  = [];
cfg.output  = 'powandcsd';
cfg.channel = 'all';
cfg.method  = 'mtmfft';
cfg.taper   = 'boxcar';
dataFreq = ft_freqanalysis(cfg, dataPre);
 
% source reconstruction
cfg             = [];
cfg.method      = 'eloreta';
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = vol.vol;
sourceFreq      = ft_sourceanalysis(cfg, dataFreq);  % compute the source model

% -----------------------------

% Compute an ERP
cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = [-1 0]; % calculate the average of the covariance matrices 
                                   % for each trial (but using the pre-event baseline  data only)
dataAvg = ft_timelockanalysis(cfg, dataPre);

% source reconstruction
cfg             = [];
cfg.method      = 'eloreta';
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = vol.vol;
sourceTime      = ft_sourceanalysis(cfg, dataAvg); 