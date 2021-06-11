function test_pull1757

%addpath('..');

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

%% get multitaper fourierspctr
cfg = [];
cfg.method     = 'mtmfft';
cfg.taper      = 'dpss';      % multitaper
cfg.output     = 'fourier';
cfg.keeptapers = 'yes';       % fourier requires keeping tapers
cfg.keeptrials = 'yes';       % and trials
cfg.pad        ='nextpow2';

cfg.foi        = 12;
cfg.tapsmofrq  = 4;
data_FFT = ft_freqanalysis(cfg, data);
% fourierspctrm size = ntaper*ntrial x nchan = 7*500 x 4


%% use new version of powcorr_ortho
cfg = [];
cfg.method = 'powcorr_ortho';
stat = ft_connectivityanalysis(cfg, data_FFT);
imagesc(stat.powcorrspctrm);


