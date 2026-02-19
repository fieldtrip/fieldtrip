function test_pull1757

% WALLTIME 00:10:00
% MEM 1gb
% DATA no

%% generate data (from test_bug2148)
cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 4;
cfg.method      = 'ar';
cfg.params(:,:,1) = [.8   0   0 .2;
    0  .9  .5  0;
    .4   0  .5 .3;
    0  .2   0 .7];
cfg.params(:,:,2) = [-.5   0   0  0;
    0 -.8   0  0;
    0   0 -.2  0;
    0   0   0 .1];
cfg.noisecov      = [.3 0  0  0;
    0 1  0  0;
    0 0 .2  0;
    0 0  0  .4];
[data] = ft_connectivitysimulation(cfg);

%% get multitaper fourierspctrm
cfg = [];
cfg.method     = 'mtmfft';
cfg.taper      = 'dpss';      % multitaper
cfg.output     = 'fourier';
cfg.keeptapers = 'yes';       % fourier requires keeping tapers
cfg.keeptrials = 'yes';       % and trials
cfg.pad        ='nextpow2';

cfg.foi        = 12;
cfg.tapsmofrq  = 4;
freq = ft_freqanalysis(cfg, data);
% fourierspctrm size = ntaper*ntrial x nchan = 7*500 x 4


%% use new version of powcorr_ortho
cfg = [];
cfg.method = 'powcorr_ortho';
stat = ft_connectivityanalysis(cfg, freq);
imagesc(stat.powcorrspctrm);

% should also work for time-resolved data
cfg = [];
cfg.method = 'mtmconvol';
cfg.toi    = (0.25:0.05:0.75);
cfg.foi    = (2:2:20);
cfg.taper  = 'hanning';
cfg.t_ftimwin = ones(1, numel(cfg.foi))./2;
cfg.output = 'fourier';
freq = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'powcorr_ortho';
stat = ft_connectivityanalysis(cfg, freq);
