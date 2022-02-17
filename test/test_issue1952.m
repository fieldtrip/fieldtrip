%% PATHS

clear variables
restoredefaultpath
addpath /home/lau/matlab/fieldtrip/
ft_defaults

path = fullfile('/home/lau/projects/functional_cerebellum/raw', ...
                '0001/20210810_000000/', ...
                'MEG/001.func_cerebellum_raw/files');

%% SEGMENT

cfg                         = [];
cfg.dataset                 = fullfile(path, 'func_cerebellum_raw.fif');
cfg.trialdef.eventtype      = 'STI101';
cfg.trialdef.eventvalue     = 3;
cfg.trialdef.prestim        = 0.200;        
cfg.trialdef.poststim       = 0.600;        
cfg.continuous              = 'yes';
cfg = ft_definetrial(cfg);

% preprocess the data
cfg.channel                 = {'MEG'};
cfg.demean                  = 'yes';     % apply baselinecorrection
cfg.baselinewindow          = [-0.200 0]; % on basis of mean signal between -0.05 and 0

data_fif = ft_preprocessing(cfg);

%% TIMELOCK

cfg = [];

timelock = ft_timelockanalysis(cfg, data_fif);

%% PLOT TIMELOCK

cfg = [];
cfg.channel = 'MEGMAG';

ft_multiplotER(cfg, timelock);

%% DENOISE WITH SSP

cfg = [];
cfg.ssp = 'all';

denoised_timelock = ft_denoise_ssp(cfg, timelock);

%% PLOT DENOISED TIMELOCK

cfg = [];
cfg.channel = 'MEGMAG';

ft_multiplotER(cfg, denoised_timelock);