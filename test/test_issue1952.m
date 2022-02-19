function test_issue1952

% WALLTIME 00:20:00
% MEM 4gb

% DEPENDENCY ft_denoise_ssp

%% PATHS

% path = '/home/lau/Dokumenter/wakeman_henson_face_data/ds117/sub001/MEG';
% dataset = fullfile(path, 'run_01_raw.fif');

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/workshop/natmeg/oddball1_mc_downsampled.fif');
if ~exist(dataset, 'file')
  % probably not on filesystem at Donders
  datadir = tempdir;
  s = ftp('ftp.fieldtriptoolbox.org');
  cd(s,'pub/fieldtrip/workshop/natmeg');
  mget(s, 'oddball1_mc_downsampled.fif', datadir);
  close(s);
  dataset = fullfile(datadir, 'oddball1_mc_downsampled.fif');
end

%% SEGMENT

cfg                         = [];
cfg.dataset                 = dataset;
cfg.trialdef.eventtype      = 'STI101';
cfg.trialdef.eventvalue     = 5;
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

ft_multiplotER(cfg, timelock);

%% DENOISE WITH SSP

cfg = [];
cfg.ssp = 'all';

denoised_timelock = ft_denoise_ssp(cfg, timelock);

%% PLOT SINGLE PLOT

cfg = [];
cfg.channel = 'MEG2611';

ft_singleplotER(cfg, denoised_timelock);
ft_singleplotER(cfg, timelock);

%% PLOT DENOISED TIMELOCK

cfg = [];

ft_multiplotER(cfg, denoised_timelock); %% fails