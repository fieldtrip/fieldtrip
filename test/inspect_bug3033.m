% function inspect_bug3033

% TEST inspect_bug3033
% TEST ft_plot_topo ft_databrowser ft_topoplotER

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3033'));

if true
  % this section only needs to run once
  cfg = [];
  cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
  cfg.continuous = 'yes';
  cfg.trl(:,1) = (1:300:10000)';
  cfg.trl(:,2) = (1:300:10000)' + 299;
  cfg.trl(:,3) = 0;
  cfg.demean = 'yes';
  ctf = ft_preprocessing(cfg);
  
  save ctf ctf
  
  cfg = [];
  cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/natmeg/oddball1_mc_downsampled.fif');
  cfg.continuous = 'yes';
  cfg.trl(:,1) = (1:1000:30000)';
  cfg.trl(:,2) = (1:1000:30000)' + 999;
  cfg.trl(:,3) = 0;
  cfg.demean = 'yes';
  elekta = ft_preprocessing(cfg);
  
  save elekta elekta
end

load ctf
load elekta

%% make some derived data from the original CTF file

cfg = [];
cfg.channel = 'MEG';
ctf_meg = ft_selectdata(cfg, ctf);

tcfg = [];
tcfg.template = 'ctf151_neighb.mat';
tcfg.method = 'template';

cfg = [];
cfg.neighbours = ft_prepare_neighbours(tcfg);
ctf_meg_planar = ft_megplanar(cfg, ctf_meg);

cfg = [];
ctf_meg_planar_timelock = ft_timelockanalysis(cfg, ctf_meg_planar);

cfg = [];
cfg.method = 'pca'; % simple and robust
ctf_meg_planar_comp = ft_componentanalysis(cfg, ctf_meg_planar);

cfg = [];
cfg.method = 'wavelet';
ctf_meg_planar_freq = ft_timelockanalysis(cfg, ctf_meg_planar);


%% make some derived data from the original Elekta file

cfg = [];
cfg.channel = 'EEG';
elekta_eeg = ft_selectdata(cfg, elekta);

cfg = [];
cfg.channel = 'MEG';
elekta_meg = ft_selectdata(cfg, elekta); % both magnetometers and planar gradiometers

cfg = [];
cfg.channel = 'megmag';
elekta_mag = ft_selectdata(cfg, elekta);

cfg = [];
cfg.channel = 'meggrad';
elekta_planar = ft_selectdata(cfg, elekta);

% start some computations

cfg = [];
elekta_eeg_timelock     = ft_timelockanalysis(cfg, elekta_eeg);
elekta_meg_timelock     = ft_timelockanalysis(cfg, elekta_meg);
elekta_mag_timelock     = ft_timelockanalysis(cfg, elekta_mag);
elekta_planar_timelock  = ft_timelockanalysis(cfg, elekta_planar);

cfg = [];
cfg.method = 'pca'; % simple and robust
elekta_eeg_comp     = ft_componentanalysis(cfg, elekta_eeg);
elekta_meg_comp     = ft_componentanalysis(cfg, elekta_meg);
elekta_mag_comp     = ft_componentanalysis(cfg, elekta_mag);
elekta_planar_comp  = ft_componentanalysis(cfg, elekta_planar);


%%


% make a copy for convenience
data = elekta_mag;

% this fails due to EEG/MEG confusion
cfg = [];
cfg.layout = ft_prepare_layout(cfg, data);

ft_databrowser(cfg, data);

