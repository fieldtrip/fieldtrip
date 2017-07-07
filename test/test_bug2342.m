function test_bug2342

% MEM 1500mb
% WALLTIME 00:10:00


% first create some data
%--------------------------------------------------------
% make 3 channels with no direct link between 1 and 2
cfg             = [];
cfg.ntrials     = 100;
cfg.triallength = 0.5;
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
data2           = data;
data2.label     = {'signal001b';'signal002b';'signal003b'};
data            = ft_appenddata([], data, data2);

% according to Martin, padding should get rid of most the zigzags.

% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 4;
freq           = ft_freqanalysis(cfgf, data);
cfgf.padding   = 10;
freqpad        = ft_freqanalysis(cfgf, data);

% connectivityanalysis
cfg = [];
cfg.method           = 'granger';
cfg.granger.sfmethod = 'bivariate';
g    = ft_connectivityanalysis(cfg, freq);
gpad = ft_connectivityanalysis(cfg, freqpad);

%...but it doesn't

% now make it more realistic and add a bit of noise to signals 4-6 (to make
% it not perfectly collinear).
data2 = data;
for k = 1:numel(data.trial)
  data2.trial{k}(4:6,:) = data.trial{k}(4:6,:)+randn(3,100)*0.000001;
end

cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 2;
freq2           = ft_freqanalysis(cfgf, data2);
cfgf.padding   = 10;
freq2pad        = ft_freqanalysis(cfgf, data2);

% connectivityanalysis
cfg = [];
cfg.method           = 'granger';
cfg.granger.sfmethod = 'bivariate';
g2    = ft_connectivityanalysis(cfg, freq2);
g2pad = ft_connectivityanalysis(cfg, freq2pad);

% conclusion: adding a tiny bit of noise removes the zigzags.
