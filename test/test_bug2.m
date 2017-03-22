function test_bug2

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqanalysis ft_megplanar 

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedaveraging/dataFC_LP.mat'));

cfg = [];
cfg.trials = 1:5;
data = ft_preprocessing(cfg, dataFC_LP);
clear dataFC_LP

cfg = [];
cfg.method = 'mtmfft';
cfg.foilim = [1 50];
cfg.taper = 'hanning';
cfg.output = 'fourier';
freq1 = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'mtmfft';
cfg.foilim = [1 50];
cfg.taper = 'dpss';
cfg.tapsmofrq = 5;
cfg.output = 'fourier';
freq2 = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'mtmconvol';
cfg.foi = 1:10;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(1,10);
cfg.toi = 0:0.1:1;
cfg.output = 'fourier';
freq3 = ft_freqanalysis(cfg, data);

% compared to the time when the bug was initially filed, ft_megplanar has been
% updated and now requires an explicit specification of cfg.neighbours
cfg = [];
cfg.method = 'template';
neighbours = ft_prepare_neighbours(cfg, data);

% convert to planar representation
cfg = [];
cfg.neighbours = neighbours;
dataP  = ft_megplanar(cfg, data);
freq1P = ft_megplanar(cfg, freq1);
freq2P = ft_megplanar(cfg, freq2);

