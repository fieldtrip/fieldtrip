function test_pull783

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/epilepsy/raw/case2/neuromag/case2_sss.fif');

% Idee, split half voor een enkel kanaaltype, sommige in m en andere in cm

%%

cfg = [];
cfg.dataset = filename;
cfg.trl = [1 1000 0];

cfg.coilaccuracy = [];
cfg.siunits = 'no';
data1 = ft_preprocessing(cfg);

cfg.coilaccuracy = [];
cfg.siunits = 'yes';
data2 = ft_preprocessing(cfg);

cfg.coilaccuracy = 1;
cfg.siunits = 'no';
data3 = ft_preprocessing(cfg);

cfg.coilaccuracy = 1;
cfg.siunits = 'yes';
data4 = ft_preprocessing(cfg);

%% do some processing

cfg = [];
timelock1 = ft_timelockanalysis(cfg, data1);
timelock2 = ft_timelockanalysis(cfg, data2);
timelock3 = ft_timelockanalysis(cfg, data3);
timelock4 = ft_timelockanalysis(cfg, data4);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'fourier';
cfg.foilim = [1 50];
freq1 = ft_freqanalysis(cfg, data1);
freq2 = ft_freqanalysis(cfg, data2);
freq3 = ft_freqanalysis(cfg, data3);
freq4 = ft_freqanalysis(cfg, data4);

%% pretend that we are estimating the noise covariance

cfg = [];
cfg.covariance = 'yes';
noise1 = ft_timelockanalysis(cfg, data1);
noise2 = ft_timelockanalysis(cfg, data2);
noise3 = ft_timelockanalysis(cfg, data3);
noise4 = ft_timelockanalysis(cfg, data4);

%% raw data

cfg = [];
cfg.split = 'all';
pw_data1 = ft_denoise_prewhiten(cfg, data1, noise1);
pw_data2 = ft_denoise_prewhiten(cfg, data2, noise2);
pw_data3 = ft_denoise_prewhiten(cfg, data3, noise3);
pw_data4 = ft_denoise_prewhiten(cfg, data4, noise4);

%% timelock data

cfg = [];
cfg.split = 'all';
pw_timelock1 = ft_denoise_prewhiten(cfg, timelock1, noise1);
pw_timelock2 = ft_denoise_prewhiten(cfg, timelock2, noise2);
pw_timelock3 = ft_denoise_prewhiten(cfg, timelock3, noise3);
pw_timelock4 = ft_denoise_prewhiten(cfg, timelock4, noise4);

%% frequency data

cfg = [];
cfg.split = 'all';
pw_freq1 = ft_denoise_prewhiten(cfg, freq1, noise1);
pw_freq2 = ft_denoise_prewhiten(cfg, freq2, noise2);
pw_freq3 = ft_denoise_prewhiten(cfg, freq3, noise3);
pw_freq4 = ft_denoise_prewhiten(cfg, freq4, noise4);


