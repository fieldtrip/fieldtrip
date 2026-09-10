function test_tutorial_denoising_opm(datadir)

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_preprocessing ft_timelockanalysis ft_multiplotER ft_singleplotER ft_topoplotER ft_denoise_hfc ft_denoise_ssp
% DATA public

% see http://www.fieldtriptoolbox.org/tutorial/preproc/denoising_opm
% this testscript corresponds to the version on the wiki at 10 September 2026

if nargin<1
  datadir = dccnpath('/project/3031000.02/test/original/meg/fieldline/sub-001/ses-opm01/meg/');
end

dataset = fullfile(datadir, 'sub-001_ses-opm01_task-emptyroom_meg.fif');

% first pass for exploration, to evaluate the power spectrum, and perform PCA
cfg         = [];
cfg.dataset = dataset;
data_er     = ft_preprocessing(cfg);

%
cfg = [];
cfg.method = 'mtmfft';
cfg.taper  = 'hanning';
cfg.channel = 'MEG';
freq_er = ft_freqanalysis(cfg, data_er);

cfg = [];
cfg.operation = 'log10';
cfg.parameter = 'powspctrm';
freq_er = ft_math(cfg, freq_er);

figure
subplot(121); plot(freq_er.freq, freq_er.powspctrm); xlim([0 80]);
subplot(122); plot(freq_er.freq, freq_er.powspctrm); xlim([0 500]);

%
cfg         = [];
cfg.length  = 2;
data_er     = ft_redefinetrial(cfg, data_er);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper  = 'hanning';
cfg.channel = 'MEG';
freq_er = ft_freqanalysis(cfg, data_er);

cfg = [];
cfg.operation = 'log10';
cfg.parameter = 'powspctrm';
freq_er = ft_math(cfg, freq_er);

figure
subplot(121); plot(freq_er.freq, freq_er.powspctrm); xlim([0 80]);
subplot(122); plot(freq_er.freq, freq_er.powspctrm); xlim([0 500]);
%
%
% The **[ft_redefinetrial](/reference/ft_redefinetrial)** function cuts the continuous data into 2-second segments, which helps reduce the memory requirements and improves the spectral estimation.
%
% We can visualize the spatial distribution of noise across sensors using a butterfly plot:
%
cfg        = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';
cfg.linecolor = 'spatial';
cfg.viewmode  = 'butterfly';
ft_multiplotER(cfg, freq_er);

%
cfg = [];
cfg.method = 'pca';
cfg.channel = 'MEG';
comp_er = ft_componentanalysis(cfg, data_er);

cfg        = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';
%[rej, art] = ft_icabrowser(cfg, comp_er);

%
cfg = [];
cfg.dataset = dataset;
cfg.bpfilter = 'yes';
cfg.bpfreq   = [1 80];
cfg.bpfilttype = 'firws';
cfg.bsfilter  = 'yes';
cfg.bsfreq    = [58 62];
cfg.bsfilttype = 'firws';
data_er = ft_preprocessing(cfg);

cfg         = [];
cfg.length  = 2;
data_er     = ft_redefinetrial(cfg, data_er);

cfg = [];
cfg.method = 'pca';
cfg.channel = 'MEG';
comp_er = ft_componentanalysis(cfg, data_er);

cfg        = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';
%[rej, art] = ft_icabrowser(cfg, comp_er);

%
dataset = fullfile(datadir, 'sub-001_ses-opm01_task-veful_meg.fif');

cfg                    = [];
cfg.dataset            = dataset;
cfg.trialdef.eventtype = 'stim_onset';
cfg.trialdef.prestim   = 0.1;
cfg.trialdef.poststim  = 0.5;
cfg = ft_definetrial(cfg);

cfg.bpfilter = 'yes';
cfg.bpfilttype = 'firws';
cfg.bpfreq = [1 80];
cfg.demean = 'yes';
cfg.padding = 4; 
cfg.baselinewindow = [-0.1 0];
cfg.bsfilter = 'yes';
cfg.bsfreq  = [58 62];
cfg.bsfilttype = 'firws';
cfg.channel = data_er.label;
data = ft_preprocessing(cfg);

%
cfg = [];
cfg.channel = 'MEG';
tlck = ft_timelockanalysis(cfg, data);

cfg        = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';
cfg.linecolor = 'spatial';
cfg.viewmode = 'butterfly';
ft_multiplotER(cfg, tlck);

%
cfg            = [];
cfg.channel    = 'MEG';
cfg.keeptrials = 'yes';
tlck_trials    = ft_timelockanalysis(cfg, data);

tlck_trials0 = tlck_trials;
tlck_trials0.trial(:) = 0;

nrpt   = numel(data.trial);
design = [ones(1,nrpt) ones(1,nrpt)*2; 1:nrpt 1:nrpt];

cfg           = [];
cfg.method    = 'analytic';
cfg.statistic = 'depsamplesT';
cfg.design    = design;
stat = ft_timelockstatistics(cfg, tlck_trials, tlck_trials0);

cfg        = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';
cfg.linecolor = 'spatial';
cfg.viewmode = 'butterfly';
cfg.parameter = 'stat';
ft_multiplotER(cfg, stat);

% compute the SSPs on-the-fly and clean the task data
cfg = [];
cfg.channel    = 'MEG';
cfg.refchannel = 'MEG';
data_ssp1 = ft_denoise_ssp(cfg, data, data);    % use the data itself for the ssp estimation
data_ssp2 = ft_denoise_ssp(cfg, data, data_er); % use the emptyroom data for the ssp estimation

% compute ERF
cfg = [];
cfg.channel = 'MEG';
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 0];
tlck_ssp1 = ft_timelockanalysis(cfg, data_ssp1);
tlck_ssp2 = ft_timelockanalysis(cfg, data_ssp2);

cfg        = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';
cfg.figure = 'subplot';
ft_topoplotER(cfg, tlck, tlck_ssp1, tlck_ssp2);

%
cfg            = [];
cfg.channel    = 'MEG';
cfg.keeptrials = 'yes';
tlck_trials_ssp1 = ft_timelockanalysis(cfg, data_ssp1);
tlck_trials_ssp2 = ft_timelockanalysis(cfg, data_ssp2);

tlck_trials0 = tlck_trials_ssp1;
tlck_trials0.trial(:) = 0;

nrpt   = numel(data.trial);
design = [ones(1,nrpt) ones(1,nrpt)*2; 1:nrpt 1:nrpt];

cfg           = [];
cfg.method    = 'analytic';
cfg.statistic = 'depsamplesT';
cfg.design    = design;
stat_ssp1 = ft_timelockstatistics(cfg, tlck_trials_ssp1, tlck_trials0);
stat_ssp2 = ft_timelockstatistics(cfg, tlck_trials_ssp2, tlck_trials0);

cfg        = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';
cfg.figure = 'subplot';
cfg.parameter = 'stat';
cfg.xlim   = [0.07 0.075];
ft_topoplotER(cfg, stat, stat_ssp1, stat_ssp2);

% clean the data with HFC
cfg = [];
cfg.order = 2;
data_hfc = ft_denoise_hfc(cfg, data);

% clean the data with AMM
cfg = [];
cfg.amm.thr = 1-1e-8;
data_amm = ft_denoise_amm(cfg, data);
