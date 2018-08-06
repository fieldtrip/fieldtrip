function test_tutorial_tmseeg

% MEM 16gb
% WALLTIME 01:20:00

% TEST ft_math ft_interpolatenan

triggers = {'S  1', 'S  3'}; % These values correspond to the markers placed in this dataset

cfg = [];
cfg.dataset                 = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/tms/sp/jimher_toolkit_demo_dataset_.eeg');
cfg.continuous              = 'yes';
cfg.trialdef.prestim        = .5;   % Data to read in prior to event onset
cfg.trialdef.poststim       = 1.5;  % Data to read in after event onset
cfg.trialdef.eventtype     = 'Stimulus' ;
cfg.trialdef.eventvalue     = triggers ;
cfg = ft_definetrial(cfg); % Create trial structure

% We can now use this trial structure (located in cfg.trl) to read our trials from disk into memory. Because we will need this trial structure later, we will save it into another variable.
trl = cfg.trl;

cfg.channel = {'all' '-5' '-mastoid L' '-mastoid R'}; % Here we indicate the channels we would like to read and/or exclude.
cfg.reref = 'yes'; % We will rereference our data
cfg.refchannel = {'all'}; % Here we specify our reference channels
cfg.implicitref = '5'; % Here we can specify the name of our implicit reference channel after rereferencing

data_tms_raw = ft_preprocessing(cfg);

%save('data_tms_raw','data_tms_raw','-v7.3');

if false
  % Inspect using databrowser, this should not run in automated mode
  cfg = [];
  cfg.preproc.demean = 'yes';
  cfg.preproc.baselinewindow = [-0.1 -0.001]; % For plotting purposes we will apply a baseline correction to the pre-stimulation period
  ft_databrowser(cfg, data_tms_raw);
end

% Average
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 -0.001];

data_tms_avg = ft_timelockanalysis(cfg, data_tms_raw);

% clear data_tms_raw to save memory
clear data_tms_raw

% plot all in seperate window
for i=1:numel(data_tms_avg.label) % Loop through all channels
  figure;
  plot(data_tms_avg.time, data_tms_avg.avg(i,:)); % Plot all data
  xlim([-0.1 0.6]); % Here we can specify the limits of what to plot on the x-axis
  ylim([-23 15]); % Here we can specify the limits of what to plot on the y-axis
  title(['Channel ' data_tms_avg.label{i}]);
  ylabel('Amplitude (uV)')
  xlabel('Time (s)');
end;

% close all windows
close all;

% Plot channel
channel = '17';

figure;
i = find(strcmp(channel, data_tms_avg.label));
plot(data_tms_avg.time, data_tms_avg.avg(i,:)); % Plot all data
xlim([-0.1 0.6]); % Here we can specify the limits of what to plot on the x-axis
ylim([-23 15]); % Here we can specify the limits of what to plot on the y-axis
title(['Channel ' data_tms_avg.label{i}]);
ylabel('Amplitude (uV)')
xlabel('Time (s)');

% adjust limits
xlim([-0 0.020]);
ylim([-60 100]);

% highlight artifacts
channel = '17';

figure;
channel_idx = find(strcmp(channel, data_tms_avg.label));
plot(data_tms_avg.time, data_tms_avg.avg(channel_idx,:)); % Plot all data
xlim([-0.1 0.6]); % Here we can specify the limits of what to plot on the x-axis
ylim([-60 100]); % Here we can specify the limits of what to plot on the y-axis
title(['Channel ' data_tms_avg.label{channel_idx}]);
ylabel('Amplitude (uV)')
xlabel('Time (s)');

%Highligh segments
hold on;

ringing = [-0.0002 0.0044];
muscle = [0.0044 0.015];
decay = [0.015 0.200];
recharge = [0.4994 0.5112];


colors = 'rgcm';
labels = {'ringing','muscle','decay','recharge'};
artifacts = [ringing; muscle; decay; recharge];

for i=1:numel(labels);
  highlight_idx = [nearest(data_tms_avg.time,artifacts(i,1)) nearest(data_tms_avg.time,artifacts(i,2)) ];
  plot(data_tms_avg.time(highlight_idx(1):highlight_idx(2)), data_tms_avg.avg(channel_idx,highlight_idx(1):highlight_idx(2)),colors(i));
end;
legend(['raw data', labels]);


%% Modified pipeline 2 - ICA on segments without ringing and recharge

% Ringing
trigger = {'S  1','S  3'};
cfg                         = [];
cfg.method                  = 'marker';
cfg.dataset                 = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/tms/sp/jimher_toolkit_demo_dataset_.eeg');
cfg.prestim                 = .001;
cfg.poststim                = .006;
cfg.trialdef.eventtype      = 'Stimulus';
cfg.trialdef.eventvalue     = trigger ;
cfg_ringing = ft_artifact_tms(cfg); % Detect TMS artifacts

% Recharge
cfg.prestim                 = -.499;
cfg.poststim                = .511;
cfg_recharge = ft_artifact_tms(cfg); % Detect TMS artifacts

% Combine into one structure
cfg_artifact = [];
cfg_artifact.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/tms/sp/jimher_toolkit_demo_dataset_.eeg');
cfg_artifact.artfctdef.ringing.artifact = cfg_ringing.artfctdef.tms.artifact; % Add ringing/step response artifact
cfg_artifact.artfctdef.recharge.artifact   = cfg_recharge.artfctdef.tms.artifact; % Add recharge artifact

% cfg_artifact = cfg_ringing;
% cfg_artifact.artfctdef.tms.artifact = [cfg_artifact.artfctdef.tms.artifact;cfg_recharge.artfctdef.tms.artifact];

cfg_artifact.artfctdef.reject        = 'partial'; % Set partial artifact rejection

cfg_artifact.trl = trl;
cfg_artifact.artfctdef.minaccepttim = 0.01;
cfg = ft_rejectartifact(cfg_artifact); % Reject trials partially

%% Read-in segmented data
cfg.channel = {'all' '-5' '-mastoid L' '-mastoid R'};
cfg.reref = 'yes';
cfg.refchannel = {'all'};
cfg.implicitref = '5';

data_tms_segmented = ft_preprocessing(cfg);


%% comp_tmsare segmented vs raw
if false
  % Inspect using databrowser, this should not run in automated mode
  
  % segmented
  cfg = [];
  cfg.artfctdef = cfg_artifact.artfctdef;
  cfg.continuous = 'yes';
  ft_databrowser(cfg, data_tms_segmented);
  
  % raw
  cfg = [];
  cfg.artfctdef = cfg_artifact.artfctdef;
  cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/tms/sp/jimher_toolkit_demo_dataset_.eeg');
  ft_databrowser(cfg);
end

close all;

%% Perform ICA on segmented data

% cfg = [];
% cfg.demean = 'yes';
% cfg.method = 'fastica';
% cfg.fastica.approach = 'symm';
% cfg.fastica.g = 'gauss';
%
% comp_tms = ft_conentanalysis(cfg, data_tms_segmented);

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/tms/sp/comp_tms.mat'));

%save('comp_tms','comp_tms','-v7.3');

%% Time-lock average comp_tmsonents
cfg = [];
cfg.vartrllength  = 2;
comp_tms_avg = ft_timelockanalysis(cfg, comp_tms);

%% Databrowser
if false
  % Inspect using databrowser, this should not run in automated mode
  figure;
  cfg = [];
  ft_databrowser(cfg, comp_tms_avg);
end

%% ft_topoplotIC
figure;
cfg = [];
cfg.component = 1:60;
cfg.comment = 'no';
cfg.layout = 'easycapM10';
ft_topoplotIC(cfg, comp_tms);

%% databrowser to browse through trials
if false
  % Inspect using databrowser, this should not run in automated mode
  cfg = [];
  cfg.layout = 'easycapM10';
  cfg.viewmode = 'comp_tmsonent';
  ft_databrowser(cfg, comp_tms);
end

%% Apply unmixing matrix to same data without demeaning

cfg = [];
cfg.demean = 'no';
cfg.unmixing = comp_tms.unmixing;
cfg.topolabel = comp_tms.topolabel;

comp_tms = ft_componentanalysis(cfg, data_tms_segmented);

% Remove comp_tmsonents
cfg = [];
cfg.component = [ 41 56 7 33 1 25 52 37 49 50 31];
cfg.demean = 'no';

data_tms_clean_segmented = ft_rejectcomponent(cfg, comp_tms);

%% Have a look at result of removal.
cfg = [];
cfg.vartrllength = 2;
cfg.preproc.demean = 'no';

data_tms_clean_avg = ft_timelockanalysis(cfg, data_tms_clean_segmented);

for i=1:numel(data_tms_clean_avg.label) % Loop through all channels
  figure;
  plot(data_tms_clean_avg.time, data_tms_clean_avg.avg(i,:),'b'); % Plot all data
  xlim([-0.1 0.6]); % Here we can specify the limits of what to plot on the x-axis
  title(['Channel ' data_tms_clean_avg.label{i}]);
  ylabel('Amplitude (uV)')
  xlabel('Time (s)');
end;


%% Restructure trials and interpolate

% Apply original structure to segmented data, gaps will be filled with nans
cfg = [];
cfg.trl = trl;
data_tms_clean = ft_redefinetrial(cfg, data_tms_clean_segmented); % Restructure cleaned data

% Replacing muscle artifact with nans
muscle_window = [0.006 0.015];
muscle_window_idx = [nearest(data_tms_clean.time{1},muscle_window(1)) nearest(data_tms_clean.time{1},muscle_window(2))];
for i=1:numel(data_tms_clean.trial)
  data_tms_clean.trial{i}(:,muscle_window_idx(1):muscle_window_idx(2))=nan;
end;

% Interpolate nans using cubic interpolation
cfg = [];
cfg.method = 'pchip';
cfg.prewindow = 0.01;
cfg.postwindow = 0.01;
data_tms_clean = ft_interpolatenan(cfg, data_tms_clean);

%% comp_tmsare
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 -0.001];

data_tms_clean_avg = ft_timelockanalysis(cfg, data_tms_clean);

for i=1:numel(data_tms_avg.label) % Loop through all channels
  figure;
  plot(data_tms_avg.time, data_tms_avg.avg(i,:),'r'); % Plot all data
  hold on;
  plot(data_tms_clean_avg.time, data_tms_clean_avg.avg(i,:),'b'); % Plot all data
  xlim([-0.1 0.6]); % Here we can specify the limits of what to plot on the x-axis
  ylim([-23 15]); % Here we can specify the limits of what to plot on the y-axis
  title(['Channel ' data_tms_avg.label{i}]);
  ylabel('Amplitude (uV)')
  xlabel('Time (s)');
  legend({'Raw' 'Cleaned'});
end;

%% Apply rest of processing steps
cfg = [];
cfg.resamplefs = 1000;
cfg.detrend = 'no';
cfg.demean = 'yes';
data_tms_clean = ft_resampledata(cfg, data_tms_clean);

%save('data_tms_clean','data_tms_clean','-v7.3');

%% Analysis - 1. time-locked average
cfg = [];
cfg.trials = find(data_tms_clean.trialinfo==1); % Find all trials corresponding to the relax condition
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.05 -0.001];
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 35;

relax_avg = ft_timelockanalysis(cfg, data_tms_clean);

cfg.trials = find(data_tms_clean.trialinfo==3);  % Find all trials corresponding to the contract condition
contract_avg = ft_timelockanalysis(cfg, data_tms_clean);

% calculate difference
cfg = [];
cfg.operation = 'subtract'; % Operation to apply
cfg.parameter = 'avg'; % The field in the data structure to which to apply the operation
difference_avg = ft_math(cfg, contract_avg, relax_avg);

% Plot TEPs of both conditions
cfg = [];
cfg.layout = 'easycapM10';
cfg.channel = '17';
cfg.xlim = [-0.1 0.6];
ft_singleplotER(cfg, relax_avg, contract_avg, difference_avg);
ylabel('Amplitude (uV)');
xlabel('time (s)');
title('Relax vs Contract');
legend({'relax' 'contract' 'contract-relax'});

%% Plotting topographies.
figure;
cfg = [];
cfg.layout = 'easycapM10';
cfg.xlim = 0:0.05:0.55;
cfg.zlim = [-2 2];
ft_topoplotER(cfg, difference_avg);

%% Analysis 3. - Global mean field power
% GMFP according to massimini is calculated aover subjects. Here we will
% calculate it within subjects. It is here defined as:
%
% GMFP(t) = sqrt(sum(Vi:k(t)-Vmean(t))^2)/K, where K is the number of
% channels. V represents the amplitude and t is time.
%  The GMFP is calculated on the time-locked average.

% Create time-locked average
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 -.001];
cfg.preproc.bpfilter = 'yes';
cfg.preproc.bpfreq = [5 100];

cfg.trials = find(data_tms_clean.trialinfo==1); % 'relax' trials
relax_avg = ft_timelockanalysis(cfg, data_tms_clean);

cfg.trials = find(data_tms_clean.trialinfo==3); % 'contract' trials
contract_avg = ft_timelockanalysis(cfg, data_tms_clean);

% GMFP calculation
cfg = [];
cfg.method = 'amplitude';
relax_gmfp = ft_globalmeanfield(cfg, relax_avg); 
contract_gmfp = ft_globalmeanfield(cfg, contract_avg); 

%Plot GMFP
figure;
plot(relax_gmfp.time, relax_gmfp.avg,'b');
hold on;
plot(contract_gmfp.time, contract_gmfp.avg,'r');
xlabel('time (s)');
ylabel('GMFP (uv^2)');
legend({'Relax' 'Contract'});
xlim([-0.1 0.6]);
ylim([0 3]); 


%% Analysis - 3. TFRs

% Calculate Induced TFRs fpor both conditions
cfg = [];
cfg.polyremoval    = 1; % Remove mean and linear trend
cfg.output         = 'pow';
cfg.method         = 'mtmconvol';
cfg.taper          = 'hanning';
cfg.foi            = 1:50;
cfg.t_ftimwin      = 0.3.*ones(1,numel(cfg.foi));
cfg.toi            = -0.5:0.05:1.5;
cfg.trials         = find(data_tms_clean.trialinfo==1);
relax_freq         = ft_freqanalysis(cfg, data_tms_clean);

cfg.trials         = find(data_tms_clean.trialinfo==3);
contract_freq      = ft_freqanalysis(cfg, data_tms_clean);

% Remove baseline
cfg = [];
cfg.baselinetype = 'relchange';
cfg.baseline = [-0.5 -0.3];
relax_freq_bc = ft_freqbaseline(cfg, relax_freq);
contract_freq_bc = ft_freqbaseline(cfg, contract_freq);

% Calculate the difference between both conditions
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
difference_freq = ft_math(cfg, contract_freq_bc, relax_freq_bc);


%% plotting
cfg = [];
cfg.xlim = [-0.1 1.0];
cfg.zlim = [-1.5 1.5];
cfg.layout = 'easycapM10';
figure;

ft_multiplotTFR(cfg, difference_freq);

%%
cfg = [];
cfg.channel = '17';
cfg.xlim = [-0.1 1.0];
cfg.zlim = [-3 3];
cfg.layout = 'easycapM10';
figure;
subplot(1,3,1);
ft_singleplotTFR(cfg, relax_freq_bc);
ylabel('Frequency (Hz)');
xlabel('time (s)');
title('Relax');

subplot(1,3,2);
ft_singleplotTFR(cfg, contract_freq_bc);
title('Contract');
ylabel('Frequency (Hz)');
xlabel('time (s)');

subplot(1,3,3);
cfg.zlim = [-1.5 1.5];
ft_singleplotTFR(cfg, difference_freq);
title('Contract - Relax');
ylabel('Frequency (Hz)');
xlabel('time (s)');

