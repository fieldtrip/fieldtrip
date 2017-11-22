function test_bug2590

% WALLTIME 00:10:00
% MEM 1500mb

% TEST ft_removetemplateartifact ft_definetrial ft_artifact_ecg ft_preprocessing ft_timelockanalysis

clear all
close all

cd(dccnpath('/home/common/matlab/fieldtrip/data'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start by constructing an averaged ECG template

cfg = [];
cfg.dataset = 'ArtifactMEG.ds';
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);  % segment data in 10 second pieces
cfg.trl = cfg.trl(1:6,:);   % only first minute of data

cfg.artfctdef.ecg.channel = {'ECG'};
cfg.artfctdef.ecg.method = 'zvalue';
cfg.artfctdef.ecg.cutoff = 1.5;
cfg.artfctdef.ecg.padding = 0.5;
% cfg.artfctdef.ecg.inspect = {'ECG'};
cfg.artfctdef.ecg.feedback = 'no';
cfg.artfctdef.ecg.pretim = 0.3;
cfg.artfctdef.ecg.psttim = 0.3;
cfg.artfctdef.ecg.mindist = 0.5;

% find the ECG segments
cfg = ft_artifact_ecg(cfg);

cfg.trl = cfg.artfctdef.ecg.artifact;
cfg.trl(:,3) = 0;
cfg.demean = 'yes';
cfg.baselinewindow = [0 0.050];

cfg.channel = 'ECG';
data_ecg = ft_preprocessing(cfg);

% here it would be wise to remove atypical trials from the data, e.g. noisy heartbeats
% cfg = [];
% data_ecg = ft_rejectvisual(cfg, data_ecg);

figure
plot(data_ecg.time{1}, data_ecg.trial{1});
title('single trial')

cfg = [];
template_ecg = ft_timelockanalysis(cfg, data_ecg);

figure
hold on
plot(template_ecg.time, template_ecg.avg)
plot(template_ecg.time, template_ecg.avg+sqrt(template_ecg.var), 'r')
plot(template_ecg.time, template_ecg.avg-sqrt(template_ecg.var), 'r')
title('average, plus/minus standard deviation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now that we have the template, we can clean the data original data

cfg = [];
cfg.artifact = data_ecg.sampleinfo;
data_ecg_clean = ft_removetemplateartifact(cfg, data_ecg, template_ecg);

figure
hold on
plot(data_ecg      .time{1}, data_ecg      .trial{1}, 'b');
plot(data_ecg_clean.time{1}, data_ecg_clean.trial{1}, 'g');
title('single trial')
legend({'original', 'cleaned'})

cfg = [];
template_ecg_clean = ft_timelockanalysis(cfg, data_ecg_clean);

figure
hold on
plot(template_ecg.time, template_ecg.avg, 'b')
plot(template_ecg_clean.time, template_ecg_clean.avg, 'g')
title('average')
legend({'original', 'cleaned'})

figure
hold on
plot(template_ecg.time, template_ecg.var, 'b')
plot(template_ecg_clean.time, template_ecg_clean.var, 'g.')
title('variance')
legend({'original', 'cleaned'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% more interesting is to clean MEG data rather than ECG data

cfg = data_ecg.cfg; % use the same ecg-triggered segmentation
cfg.channel = {'MEG', 'ECG'};
data_meg = ft_preprocessing(cfg);

cfg = [];
template_meg = ft_timelockanalysis(cfg, data_meg);

cfg = [];
cfg.artifact = data_ecg.sampleinfo;
data_meg_clean = ft_removetemplateartifact(cfg, data_meg, template_meg);

figure
hold on
plot(data_meg      .time{1}, data_meg      .trial{1}(1,:), 'b');
plot(data_meg_clean.time{1}, data_meg_clean.trial{1}(2,:), 'g');
title('single trial')
legend({'original', 'cleaned'})

cfg = [];
cfg.channel = 'MEG';
timelock_meg       = ft_timelockanalysis(cfg, data_meg);
timelock_meg_clean = ft_timelockanalysis(cfg, data_meg_clean);

figure
hold on
plot(timelock_meg.time, timelock_meg.avg, 'b')
plot(timelock_meg_clean.time, timelock_meg_clean.avg, 'g.')
title('average')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and even more interesting is to clean MEG data that is segmented completely differently

cfg = [];
cfg.dataset = 'ArtifactMEG.ds';
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);  % segment data in 10 second pieces
cfg.trl = cfg.trl(1:6,:);   % only first minute of data
cfg.channel = {'MEG', 'ECG'};
continuous_meg = ft_preprocessing(cfg);

cfg = [];
cfg.artifact = data_ecg.sampleinfo;
continuous_meg_clean = ft_removetemplateartifact(cfg, continuous_meg, template_meg);

% there is a change in the channel ordering, which is annoying for plotting
ecgchan1 = find(strcmp(continuous_meg.label, 'ECG'));
ecgchan2 = find(strcmp(continuous_meg_clean.label, 'ECG'));

figure
hold on
plot(continuous_meg      .time{1}, continuous_meg      .trial{1}(ecgchan1,:), 'b');
plot(continuous_meg_clean.time{1}, continuous_meg_clean.trial{1}(ecgchan2,:), 'g');
title('single trial')
legend({'original', 'cleaned'})

warning('note that the first artifact was not detected by ft_artifact_ecg and hence not cleaned')
