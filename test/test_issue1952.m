function test_issue1952

% WALLTIME 00:20:00
% MEM 4gb

% DEPENDENCY ft_denoise_ssp

%% PATHS

% path = '/home/lau/Dokumenter/wakeman_henson_face_data/ds117/sub001/MEG';
% dataset = fullfile(path, 'run_01_raw.fif');

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/epilepsy/raw/case1/neuromag/case1_cHPI_raw.fif');
if ~exist(dataset, 'file')
  % probably not on filesystem at Donders
  datadir = tempdir;
  s = ftp('ftp.fieldtriptoolbox.org');
  cd(s,'pub/fieldtrip/tutorial/epilepsy/raw/case1/neuromag');
  mget(s, 'case1_cHPI_raw.fif', datadir);
  close(s);
  dataset = fullfile(datadir, 'case1_cHPI_raw.fif');
end

%% SEGMENT

cfg                         = [];
cfg.dataset                 = dataset;

hdr = ft_read_header(dataset, 'checkmaxfilter', false);
fs  = hdr.Fs;

cfg.trl = [(fs:fs:100*fs)+1;(2*fs:fs:101*fs)]';
cfg.trl(:,3) = 0;


% cfg.trialdef.eventtype      = 'STI101';
% cfg.trialdef.eventvalue     = 5;
% cfg.trialdef.prestim        = 0.200;        
% cfg.trialdef.poststim       = 0.600;        
% cfg.continuous              = 'yes';
% cfg = ft_definetrial(cfg);

% preprocess the data
cfg.channel                 = {'MEG'};
cfg.demean                  = 'yes';     % apply baselinecorrection
%cfg.baselinewindow          = [-0.200 0]; % on basis of mean signal between -0.05 and 0
cfg.checkmaxfilter = false;
data_fif = ft_preprocessing(cfg);

%% TIMELOCK

cfg = [];

timelock = ft_timelockanalysis(cfg, data_fif);

%% PLOT TIMELOCK

cfg = [];
cfg.layout = 'neuromag306mag_helmet.mat';
ft_multiplotER(cfg, timelock);

%% DENOISE WITH SSP

cfg = [];
cfg.ssp = 'all';

denoised_timelock = ft_denoise_ssp(cfg, timelock);

%% DENOISE WITH SSP CELL ARRAY

cfg = [];
% cfg.ssp = {'mag_ssp_upright_fif___pca_mags_v1', ...
%            'mag_ssp_upright_fif___pca_mags_v2', ...
%            'mag_ssp_upright_fif___pca_mags_v3'};
cfg.ssp = {'magssp68iason_fif___pca_v1', ...
  'magssp68iason_fif___pca_v2'};

denoised_timelock_specific = ft_denoise_ssp(cfg, timelock);

%% PLOT SINGLE PLOT

cfg = [];
cfg.channel = 'MEG2611';
cfg.layout = 'neuromag306mag_helmet.mat';

ft_singleplotER(cfg, timelock, denoised_timelock, denoised_timelock_specific);

%% PLOT DENOISED TIMELOCK

cfg = [];
cfg.layout = 'neuromag306mag_helmet.mat';
ft_multiplotER(cfg, timelock, denoised_timelock, denoised_timelock_specific); %% fails
