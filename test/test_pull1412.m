function test_pull1412

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_heartrate

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/pull1412'));

%%
% this corresponds to the preprocessed dataset 006_3013065.02_rest1 from bug3433

load datappg

cfg = [];
cfg.channel = 'HR';
cfg.threshold = 0.7;
cfg.method = 'findpeaks';
% cfg.method = 'pantompkin'; this does not work very well on the PPG data
heartrate0 = ft_heartrate(cfg, data);

figure
plot(heartrate0.time{1}, heartrate0.trial{1}(1,:), '-');

%%
% this corresponds to the ECG channel from ArtifactMEG.ds as documented on
% http://www.fieldtriptoolbox.org/example/use_independent_component_analysis_ica_to_remove_ecg_artifacts/

load dataecg

cfg = [];
cfg.channel = 'ECG';
cfg.threshold = 1.2;
cfg.method = 'findpeaks';
cfg.flipsignal = 'no';
heartrate1 = ft_heartrate(cfg, data);

cfg = [];
cfg.channel = 'ECG';
cfg.method = 'pantompkin';
heartrate2 = ft_heartrate(cfg, data);

figure
plot(heartrate1.time{1}, heartrate1.trial{1}(1,:), 'b-');
hold on
plot(heartrate2.time{1}, heartrate2.trial{1}(1,:), 'rx');
legend({'findpeaks', 'pantompkin'});