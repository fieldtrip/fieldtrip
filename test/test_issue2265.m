function test_issue2265

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_freqanalysis ft_specest_mtmfft
% DATA no

data = [];
data.trial{1} = randn(1,481);
data.time{1}  = linspace(-0.2, 0.6, 481);
data.label    = {'chan01'};

cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper  = 'hanning';
cfg.pad    = 1;
freq       = ft_freqanalysis(cfg, data);

% the original issue is that due to numeric rounding error the estimated
% number of padding samples was one too many, causing a mismatch between
% the actual returned frequency bins, and the specified bins in the
% freq.freq axis. This could be evaluated by observing a complex-valued
% fourier coefficient at the Nyquist frequency bin.
assert(abs(imag(freq.fourierspctrm(end)))<100*eps);
