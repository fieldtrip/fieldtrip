function test_pull1456

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_sourceanalysis ft_eloreta

load(dccnpath('/home/common/matlab/fieldtrip/data/test/test_pull1456.mat'));

% Compute spectrum
cfg         = [];
cfg.output  = 'powandcsd';
cfg.channel = 'all';
cfg.method  = 'mtmfft';
cfg.taper   = 'hanning';%'boxcar';
dataFreq1   = ft_freqanalysis(cfg, dataPre);
cfg.output  = 'fourier';
dataFreq2   = ft_freqanalysis(cfg, dataPre);

% source reconstruction
cfg             = [];
cfg.method      = 'eloreta';
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = vol.vol;
% source_eloreta1     = ft_sourceanalysis(cfg, data_eloreta1);  % compute the source model
% assert(isequal(size(source_eloreta1.avg.pow),[size(source_eloreta1.pos,1) numel(source_eloreta1.freq)]));
% 
% source_eloreta2     = ft_sourceanalysis(cfg, data_eloreta2);  % compute the source model
% i1 = find(source_eloreta2.inside,1,'first');
% assert(isequal(size(source_eloreta2.avg.mom{i1}),[3 80 193]));
% 
cfg.snr    = 10;
% cfg.method = 'mne';
% source_mne1 = ft_sourceanalysis(cfg, dataFreq1);
% source_mne2 = ft_sourceanalysis(cfg, dataFreq2);

% -> this does not work. AT ALLedit f
%cfg.method = 'harmony';
%source_harmony1 = ft_sourceanalysis(cfg, dataFreq1);
%source_harmony2 = ft_sourceanalysis(cfg, dataFreq2);

cfg.method = 'music';
source_music1 = ft_sourceanalysis(cfg, dataFreq1);
source_music2 = ft_sourceanalysis(cfg, dataFreq2);

cfg.method = 'rv';
source_rv1 = ft_sourceanalysis(cfg, dataFreq1);
source_rv2 = ft_sourceanalysis(cfg, dataFreq2);

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
