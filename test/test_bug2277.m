function test_bug2277

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preprocessing ft_preproc_padding preproc

ntime = 50;
nfreq = 30;


time_idx = unique(round(ntime/3:ntime/2));
freq_idx = unique(round(nfreq/3:nfreq/2));
freq = [];
freq.dimord = 'chan_freq_time';
freq.freq = linspace(1, 32, nfreq);
freq.time = linspace(0, 1, ntime);
freq.label = {'1'};
nchans = numel(freq.label);
freq.powspctrm = randn(nchans,nfreq, ntime);
freq.powspctrm(1, freq_idx, time_idx) = randn(1, numel(freq_idx), numel(time_idx))+4;

freq.mask = false(size(freq.powspctrm));
freq.mask(1, freq_idx, time_idx) = true;

% no masking
%cfg.maskstyle      = style used to masking, 'opacity' or 'saturation' (default = 'opacity')
%                         use 'saturation' when saving to vector-format (like *.eps) to avoid all sorts of image-problems

cfg = [];
cfg.zlim = 'maxabs';
figure;
ft_singleplotTFR(cfg, freq);

% masking options
cfg.maskparameter = 'mask';
cfg.maskalpha = 0.5;

cfg.maskstyle = 'opacity';
figure;
ft_singleplotTFR(cfg, freq);

cfg.maskstyle = 'saturation';
figure;
ft_singleplotTFR(cfg, freq);

% new contour style
cfg.maskstyle = 'outline';
figure;
ft_singleplotTFR(cfg, freq);

