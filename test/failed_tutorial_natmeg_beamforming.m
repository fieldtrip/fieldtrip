function failed_tutorial_natmeg_beamforming

% WALLTIME 00:20:00
% MEM 8gb

% this script executes the MATLAB content from
% http://www.fieldtriptoolbox.org/tutorial/natmeg/beamforming
%
% it corresponds to the wiki version of 7 October 2014

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

clear all
close all

cd(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/natmeg'));

load timefrequency/data_clean_MEG_responselocked.mat;
load dipolefitting/headmodel_meg.mat;

% Select time window of interest
cfg = [];
cfg.toilim = [0.35 0.85];
data_timewindow = ft_redefinetrial(cfg,data_clean_MEG_responselocked);

% Freqanalysis for beamformer
cfg = [];
cfg.channel      = {'MEG*2', 'MEG*3'};
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.output       = 'powandcsd';
cfg.keeptrials   = 'no';
cfg.foi          = 18;
cfg.tapsmofrq    = 4;

% for common filter over conditions
powcsd_all      = ft_freqanalysis(cfg, data_timewindow);

% for conditions
cfg.trials       = find(data_timewindow.trialinfo(:,1) == 256);
powcsd_left      = ft_freqanalysis(cfg, data_timewindow);
cfg.trials       = find(data_timewindow.trialinfo(:,1) == 4096);
powcsd_right     = ft_freqanalysis(cfg, data_timewindow);

% Create leadfield grid
cfg                 = [];
cfg.channel         = {'MEG*2', 'MEG*3'};
cfg.grad            = powcsd_all.grad;
cfg.vol             = headmodel_meg;
cfg.dics.reducerank = 2; % default for MEG is 2, for EEG is 3
cfg.grid.resolution = 0.5;   % use a 3-D grid with a 0.5 cm resolution
cfg.grid.unit       = 'cm';
cfg.grid.tight      = 'yes';
[grid] = ft_prepare_leadfield(cfg);

cfg              = [];
cfg.channel      = {'MEG*2', 'MEG*3'};
cfg.method       = 'dics';
cfg.frequency    = 18;
cfg.grid         = grid;
cfg.vol          = headmodel_meg;
cfg.senstype     = 'MEG'; % Must me 'MEG', although we only kept MEG channels, information on EEG channels is still present in data
cfg.dics.keepfilter   = 'yes'; % We wish to use the calculated filter later on
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda  = '5%';
source_all = ft_sourceanalysis(cfg, powcsd_all);

cfg              = [];
cfg.channel      = {'MEG*2', 'MEG*3'};
cfg.method       = 'dics';
cfg.frequency    = 18;
cfg.grid         = grid;
cfg.grid.filter  = source_all.avg.filter;
cfg.vol          = headmodel_meg;
cfg.senstype     ='MEG';

source_left = ft_sourceanalysis(cfg, powcsd_left);
source_right = ft_sourceanalysis(cfg, powcsd_right);

load dipolefitting/mri_realigned2.mat;

mri_resliced = ft_volumereslice([], mri_realigned2);

cfg            = [];
cfg.parameter = 'avg.pow';
source_left_int  = ft_sourceinterpolate(cfg, source_left, mri_resliced);
source_right_int  = ft_sourceinterpolate(cfg, source_right, mri_resliced);

source_diff_int  = source_left_int;
source_diff_int.avg.pow  = (source_left_int.avg.pow - source_right_int.avg.pow) ./ (source_left_int.avg.pow + source_right_int.avg.pow);

% plot
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'avg.pow';
cfg.funcolorlim   = 'maxabs';
cfg.opacitylim    = [0 1e-4];
cfg.opacitymap    = 'rampup';

ft_sourceplot(cfg, source_left_int);

cfg.location = [35 -13 76];
ft_sourceplot(cfg, source_diff_int);

load dipolefitting/headmodel_eeg.mat
load timefrequency/data_clean_EEG_responselocked.mat

% select time window
cfg = [];
cfg.toilim = [0.35 0.85];
data_timewindow = ft_redefinetrial(cfg,data_clean_EEG_responselocked);

% Freqanalysis for beamformer
cfg = [];
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.output       = 'powandcsd';
cfg.keeptrials   = 'no';
cfg.foi          = 18;
cfg.tapsmofrq    = 4;

% for common filter over conditions and full duration
powcsd_all      = ft_freqanalysis(cfg, data_timewindow);

% for conditions
cfg.trials       = find(data_timewindow.trialinfo(:,1) == 256);
powcsd_left      = ft_freqanalysis(cfg, data_timewindow);
cfg.trials       = find(data_timewindow.trialinfo(:,1) == 4096);
powcsd_right     = ft_freqanalysis(cfg, data_timewindow);

% common grid/filter
cfg                 = [];
cfg.elec            = powcsd_all.elec;
cfg.vol             = headmodel_eeg;
cfg.reducerank      = 3; % default is 3 for EEG, 2 for MEG
cfg.grid.resolution = 0.5;   % use a 3-D grid with a 0.5 cm resolution
cfg.grid.unit       = 'cm';
cfg.grid.tight      = 'yes';
[grid] = ft_prepare_leadfield(cfg);

% beamform common filter
cfg              = [];
cfg.method       = 'dics';
cfg.frequency    = 18;
cfg.grid         = grid;
cfg.vol          = headmodel_eeg;
cfg.senstype     = 'EEG'; % Remember this must be specified as either EEG, or MEG
cfg.dics.keepfilter   = 'yes';
cfg.dics.lambda       = '15%';
source_all = ft_sourceanalysis(cfg, powcsd_all);

% beamform conditions
cfg              = [];
cfg.method       = 'dics';
cfg.frequency    = 18;
cfg.grid         = grid;
cfg.grid.filter  = source_all.avg.filter; % Use the common filter
cfg.vol          = headmodel_eeg;
cfg.senstype     = 'EEG';

source_left = ft_sourceanalysis(cfg, powcsd_left);
source_right = ft_sourceanalysis(cfg, powcsd_right);

cfg            = [];
cfg.parameter = 'avg.pow';
source_left_int  = ft_sourceinterpolate(cfg, source_left, mri_resliced);
source_right_int  = ft_sourceinterpolate(cfg, source_right, mri_resliced);

source_diff_int  = source_left_int;
source_diff_int.avg.pow  = (source_left_int.avg.pow - source_right_int.avg.pow) ./ (source_left_int.avg.pow + source_right_int.avg.pow);

cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'avg.pow';
cfg.funcolorlim   = 'maxabs';

ft_sourceplot(cfg, source_left_int);

cfg.location = [-19.5 -18.5 70.5];
ft_sourceplot(cfg, source_diff_int);
