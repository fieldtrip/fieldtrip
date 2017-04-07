function test_bug1425

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_connectivityanalysis
% TEST ft_connectivity_corr

% the bug pertains to a non-specific error when trying to do coherence computation on single trial data

% reproduce the issue
data = [];
data.trial{1} = randn(3,100);
data.time{1}  = (0:99)./100;
data.label    = {'chan01';'chan02';'chan03'};

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper  = 'hanning';
freq = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'coh';
cfg.channelcmb = {'chan01' 'chan02'};
coh = ft_connectivityanalysis(cfg, freq);

% issue has been reproduced. It can be tracked down to ft_checkdata (fixcsd) where the singleton first dimension is not removed, leading to an incorrect dimensionality of the numeric data with respect to the dimord
