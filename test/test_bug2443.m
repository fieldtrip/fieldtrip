function test_bug2443

% WALLTIME 00:20:00
% MEM 4gb

% TEST ft_multiplotER

% get some data
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_neuromag306.mat');
load(filename);

data_pre  = data;
data_post = data;

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'powandcsd'; % the CSD is needed for source reconstruction
cfg.taper  = 'dpss';
cfg.foi    = 5:5:100;
cfg.tapsmofrq = 10;       % we apply plenty of frequency smoothing
freq_pre  = ft_freqanalysis(cfg, data_pre);
freq_post = ft_freqanalysis(cfg, data_post);

cfg = [];
cfg.layout      = 'neuromag306planar.lay';
cfg.parameter   = 'powspctrm';
ft_multiplotER(cfg, freq_pre, freq_post);

