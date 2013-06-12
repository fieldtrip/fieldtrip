function test_bug2148

% TEST test_bug2148
% TEST ft_connectivitysimulation ft_freqanalysis ft_connectivityanalysis
% TEST ft_connectivityplot


cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';
cfg.params(:,:,1) = [ 0.8 0   0; 
                      0   0.9 0.5;
                      0.4 0   0.5];
cfg.params(:,:,2) = [-0.5    0  0; 
                        0 -0.8  0; 
                        0    0 -0.2];
cfg.noisecov      = [0.3 0 0;
                       0 1 0;
                       0 0 0.2];

data            = ft_connectivitysimulation(cfg);

% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 2;
freq           = ft_freqanalysis(cfgf, data);

% connectivityanalysis
cfgc           = [];
cfgc.channelcmb    = {freq.label{1} freq.label{2}; freq.label{1} freq.label{3}};
cfgc.method    = 'coh';
c1             = ft_connectivityanalysis(cfgc, freq);

ft_connectivityplot([], c1)